#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <future>
#include <iostream>

#ifdef WIN32
      // For GetTickCount
    #include <Windows.h>
#endif

struct SGene
{
    enum EType { POSITIVE, NEGATIVE };

    EType type = POSITIVE;
	char name[3];
    unsigned int mask = 0;

    bool IsSame( const SGene& other ) const
    {
        return name[0] == other.name[0] && name[1] == other.name[1];
    }

    bool IsExactlySame(const SGene& other) const
    {
        return type == other.type && IsSame(other);
    }
};

struct SMutagen
{
	std::string name;

private:
	std::vector<SGene> genes;

public:
    const std::vector<SGene>& GetGenes() const { return genes; }

    void AddGene( const SGene& g )
    {
        genes.push_back(g);
        if ( g.type == SGene::NEGATIVE )
            negativeGenes |= g.mask;
        else
            positiveGenes |= g.mask;
    }

    void RemoveGene( const SGene& g )
    {
        auto iter = std::find_if(genes.begin(), genes.end(), [&g](const SGene& gene) { return g.IsSame(gene); });
        genes.erase(iter);

        negativeGenes &= ~g.mask;
        positiveGenes &= ~g.mask;
    }

    void ClearGenes()
    {
        genes.clear();
        negativeGenes = 0;
        positiveGenes = 0;
    }

    void RemoveNegativeGenes()
    {
        genes.erase( std::remove_if( genes.begin(), genes.end(), [](const SGene& g) { return g.type == SGene::NEGATIVE; } ), genes.end() );
    }

    unsigned int positiveGenes = 0;
    unsigned int negativeGenes = 0;


    void PrintGenes()
    {
        for (int i = 0; i < (int)genes.size(); ++i)
        {
            if (genes[i].type == SGene::NEGATIVE)
                printf("-");
            printf("%s", genes[i].name);
            printf(" ");
        }
        printf("\n");
    }
};

struct SSolution
{
    std::vector<const SMutagen*> steps;
    SMutagen result;
    int fitness = 0;

    std::vector<SGene> toDelete;
    std::vector<SGene> toAdd;

    SSolution()
    {
        toDelete.reserve(10);
        toAdd.reserve(10);
    }

    bool AddMutagen(const SMutagen* m)
    {
        if ( !steps.empty() && m == steps.back() )
            return false;

        toDelete.clear();
        toAdd.clear();

        for( const auto& gene : m->GetGenes() )
        {
            if ( gene.type == SGene::POSITIVE && ( result.negativeGenes & gene.mask ) )
            {
                result.RemoveGene(gene);
            }
            else if ( gene.type == SGene::NEGATIVE && ( result.positiveGenes & gene.mask ) )
            {
                result.RemoveGene(gene);
            }
            else if ( !( gene.mask & result.positiveGenes || gene.mask & result.negativeGenes ) )
            {
                result.AddGene(gene);
            }
        }
        
        steps.push_back(m);
        return true;
    }


    void PrintSteps()
    {
        for (int i = 0; i < (int)steps.size(); ++i)
            printf("%s\n", steps[i]->name.c_str());
    }
};

class Solver
{
public:
    bool SetGoal( const std::string& genes )
    {
        m_goal.name = "Exitus-1";
        if (!ParseMutagen(m_goal, genes))
            return false;

        printf("Set goal: "); m_goal.PrintGenes();

        return true;
    }    

    bool AddMutagen( const std::string& name, const std::string& genes )
    {
        SMutagen m;
        m.name = name;
        if ( !ParseMutagen(m, genes) )
            return false;

        printf("Added mutagen %s: ", name.c_str()); m.PrintGenes();

        m_mutagens.push_back(m);
        return true;
    }

    void SetThreadsCount( int threadsCount )
    {
        m_threadsCount = threadsCount;
    }

    bool Solve( SSolution& solution )
    {
        for ( int stepsCount = 1; stepsCount < 10; ++stepsCount )
        {
            printf("Seeking %i-step solutions\n", stepsCount);
            const int variantsCount = (int)pow(m_mutagens.size(), stepsCount);

            int chunk = variantsCount / m_threadsCount;
            if ( chunk < 1)
                chunk = 1;

            int start = 0;
            bool stopAll = false;
            std::vector<std::future<SSolution>> futures;
            do
            {
                int end = start + chunk;
                futures.push_back( std::async( [start, end, stepsCount, this, &stopAll]()
                        {
                            SSolution s;
                            s.steps.reserve(stepsCount);

                            for (int v = start; v < end; ++v)
                            {
                                if ( stopAll )
                                    return SSolution();

                                s.steps.clear();
                                s.result.ClearGenes();

                                bool badSeq = false;
                                for (int step = 0; step < stepsCount; ++step)
                                {
                                    const int index = (int)(v / pow(m_mutagens.size(), step)) % m_mutagens.size();
                                    if ( !s.AddMutagen(&m_mutagens[index]) )
                                    {
                                        badSeq = true;
                                        break;
                                    }
                                }

                                if ( badSeq )
                                    continue;

                                s.result.RemoveNegativeGenes();

                                if (s.result.GetGenes().size() != m_goal.GetGenes().size())
                                    continue;

                                bool allOK = true;

                                int i = 0, j = 0;
                                for ( int i = 0; i < (int)m_goal.GetGenes().size(); ++i )
                                {
                                    if (!m_goal.GetGenes()[i].IsExactlySame(s.result.GetGenes()[i]))
                                    {
                                        allOK = false;
                                        break;
                                    }
                                }

                                if (!allOK)
                                    continue;

                                s.result.PrintGenes();
                                s.PrintSteps();

                                return s;
                            }

                            return SSolution();
                        }
                    ));
                
                start = end;
            }while( start < variantsCount );
            
            while(true)
            {
                bool allReady = true;            
                for (auto& f : futures)
                {
                    if ( !f.valid() )
                        continue;

                    std::future_status status = f.wait_for(std::chrono::seconds(1));
                    if (status != std::future_status::ready)
                        allReady = false;
                    else
                    {
                        SSolution s = f.get();
                        if (s.steps.size() > 0)
                        {
                            stopAll = true;
                            return true;
                        }
                    }
                }

                if (allReady)
                    break;
            }
        }

        return false;
    }

    bool InputFromFile( const char* filename )
    {
        FILE *f = fopen(filename, "r");
        if (!f)
        {
            printf("Unable to open file %s\n", filename);
            return false;
        }

        int line = 1;
        while(!feof(f))
        {
            char mutagenName[256];
            char mutagenString[256];

            if ( fscanf(f, "%[^:]:%*[ ]%[^\n]\n", mutagenName, mutagenString) != 2 )
            {
                printf("Wrong file format at line %i\n", line );
                return false;
            }            

            if ( line == 1 )
            {
                if ( !SetGoal(mutagenString))
                    return false;
            }
            else
            {
                if ( !AddMutagen(mutagenName, mutagenString) )
                    return false;
            }

            ++line;
        }

        return true;
    }

    bool InputFromStdin()
    {
        printf("\nPlease input each know mutagen when prompted.\n");
        printf("For unknown/unwanted mutagens, just press enter\n");
        
        InputMutagen( "Exitus-1 (goal compound)", true );

        InputMutagen( "Ovid-1" );
        InputMutagen( "Ovid-2" );
        InputMutagen( "Ovid-3" );

        InputMutagen( "Solis-1" );
        InputMutagen( "Solis-2" );

        InputMutagen( "Echo-1" );
        InputMutagen( "Echo-2" );
        InputMutagen( "Echo-3" );
        InputMutagen( "Echo-4" );

        InputMutagen( "Io-1" );
        InputMutagen( "Io-2" );
        InputMutagen( "Io-3" );

        InputMutagen( "Helicon-1" );
        InputMutagen( "Helicon-2" );
        InputMutagen( "Helicon-3" );

        return true;
    }

private:
    bool InputMutagen( const char* name, bool isGoal = false )
    {
        printf("%s: ", name);
        std::string mutagenString;
        std::getline(std::cin, mutagenString);
        if ( mutagenString.empty() )
            return true;

        if ( isGoal )
            return SetGoal(mutagenString);
        else
            return AddMutagen(name, mutagenString);
    }

private:
    bool ParseMutagen( SMutagen &m, const std::string& genes )
    {        
        size_t prevPos = 0;
        size_t pos = 0;
        while ((pos = genes.find(' ', pos)) != std::string::npos)
        {
            SGene gene;
            if (!ParseGene(gene, genes, prevPos, pos))
            {
                fprintf(stderr, "Failed to parse gene at %i: %s\n", prevPos, &genes[prevPos]);
                return false;
            }

            m.AddGene(gene);
            ++pos;
            prevPos = pos;
        }

        SGene gene;
        if (!ParseGene(gene, genes, prevPos, genes.size() - 1))
        {
            fprintf(stderr, "Failed to parse gene at %i: %s\n", prevPos, &genes[prevPos]);
            return false;
        }
        m.AddGene(gene);

        return true;
    }

    bool ParseGene( SGene& gene, const std::string& genes, size_t from, size_t to )
    {
        size_t pos = from;
        if ( genes[pos] == ' ' ) ++pos;
        if ( pos > to )
        {
            fprintf(stderr, "Gene sequence is too short. Should be [-]AB\n");
            return false;
        }

        if ( genes[pos] == '-' )
        {
            gene.type = SGene::NEGATIVE;
            ++pos;
        }
        if (pos > to)
        {
            fprintf(stderr, "Gene sequence is too short. Should be [-]AB\n");
            return false;
        }

        gene.name[0] = genes[pos++];
        if (pos > to)
        {
            fprintf(stderr, "Gene sequence is too short. Should be [-]AB\n");
            return false;
        }

        gene.name[1] = genes[pos];
        gene.name[2] = '\0';

        bool known = false;
        for ( int i = 0; i < (int)m_knownGenes.size(); ++i )
        {
            if ( m_knownGenes[i].name[0] == gene.name[0] && m_knownGenes[i].name[1] == gene.name[1] )
            {
                gene.mask = m_knownGenes[i].mask;
                known = true;
                break;
            }
        }

        if ( !known )
        {
            SGeneInfo newGene;
            newGene.name[0] = gene.name[0];
            newGene.name[1] = gene.name[1];
            newGene.name[2] = 0;
            newGene.mask = 1 << m_knownGenes.size();
            m_knownGenes.push_back(newGene);
            
            gene.mask = newGene.mask;
        }

        return true;
    }

private:
    std::vector<SMutagen> m_mutagens;
    SMutagen m_goal;

    struct SGeneInfo { char name[3]; unsigned int mask; };
    std::vector<SGeneInfo> m_knownGenes;

    int m_threadsCount;
};

int main( int argc, char** argv )
{
	printf("\nUnderrail Mutagen Puzzle Solver v1.0 by MaxEd http://zxstudio.org\n");

    if ( argc < 2 )
    {
        printf("\nUsage:\n");
        printf("\n");
        printf("-f FILE -- read mutagens from file\n");
        printf("-i      -- input mutagens by hand\n");
        printf("-t N    -- use N threads (4 by default)\n");
        printf("\n");
        printf("If -f is specified, the file should have this format:\n");
        printf("\n");
        printf("MutagenName1: AA BB -CC DD EE -FF\n");
        printf("MutagenName2: AA BB -CC DD EE -FF\n");
        printf("\n");
        printf("The first mutagen in file is considered to be the goal.\n");
        printf("\n");
        printf("Example: \n");
        printf("\n");
        printf("Exitus-1: WU JJ RJ LX RU IB LM RA D2 LS CI I5 DL IQ OY\n");
        printf("Ovid-1: LX CW WU -RJ\n");
        printf("Echo-2: P9 CI OY LS OC DL RJ -CW -IQ -WU\n");
        printf("etc... (see example.txt for a complete example)\n");
        printf("\n");
        printf("Hint: to make search go faster, filter out impossible mutagens by hand,\n");
        printf("for example those containing genes that are not in the final compound\n");
        printf("and do not have a negative gene in any other sequence.\n");
        return -1;
    }

    bool interactive = false;
    std::string filename;
    int numThreads = 4;

    for( int i = 1; i < argc; ++i )
    {
        std::string arg = argv[i];
        if ( arg[0] != '-')
        {
            printf("Bad command-line argument: %s\n", arg.c_str());
            return -1;
        }

        if ( arg == "-f" )
        {
            if ( i == argc - 1)
            {
                printf("-f requires filename\n");
                return -1;
            }

            filename = argv[i+1];
            ++i;
            continue;
        }

        if ( arg == "-i" )
        {
            interactive = true;
            continue;
        }

        if ( arg == "-t" )
        {
            if (i == argc - 1)
            {
                printf("-t requires number of threads\n");
                return -1;
            }

            numThreads = atoi(argv[i+1]);

            if ( numThreads < 1 )
            {
                printf("-t requires POSITIVE number of threads\n");
                return -1;
            }

            if (numThreads > 16)
            {
                printf("Unusually large number of threads specified (%i). I assume you know what you are doing\n", numThreads);
            }

            ++i;
            continue;
        }
    }

    if ( interactive && !filename.empty() )
    {
        printf("Both -f and -i are specified. Make up your mind, dammit!\n");
        return -1;
    }

	Solver solver;
    solver.SetThreadsCount(numThreads);

    if ( !filename.empty() )
    {
        if ( !solver.InputFromFile(filename.c_str()) )
            return -1;
    }
    else
    {
        if ( !solver.InputFromStdin() )
            return -1;
    }

    //solver.AddMutagen("Ovid-1", "LX CW WU -RJ");
    ////solver.AddMutagen("Ovid-2", "JJ P9 IQ LX RU -D2");
    //solver.AddMutagen("Ovid-3", "RJ CN IQ OC -RA -CI -LX");


    ////solver.AddMutagen("Solis-1", "WU ED RA IB -JJ -I5");
    //solver.AddMutagen("Solis-2", "D2 LS CI I5 DL -OY -CW");

    //solver.AddMutagen("Echo-1", "RU IB D2");
    ////solver.AddMutagen("Echo-2", "P9 CI OY LS OC DL RJ -CW -IQ -WU");
    //solver.AddMutagen("Echo-3", "CN RU IB");
    //solver.AddMutagen("Echo-4", "IB D2 WU OY CW -RJ");

    //solver.AddMutagen("Io-1", "I5 IB IQ OY");
    ////solver.AddMutagen("Io-2", "P9 CN I5 OC -LS");
    //solver.AddMutagen("Io-3", "RJ LM RA -OC -CN -D2");

    //solver.AddMutagen("Helicon-1", "JJ D2 OY RJ LX CW -OC -IB");
    ////solver.AddMutagen("Helicon-2", "ED JJ OY I5 P9 -LM -WU -CI");
    //solver.AddMutagen("Helicon-3", "CW D2 WU -IQ");

    //solver.SetGoal("WU JJ RJ LX RU IB LM RA D2 LS CI I5 DL IQ OY");

#ifdef WIN32
    int t1 = GetTickCount();
#endif
    SSolution solution;
    if ( solver.Solve( solution ) )
    {
        printf("\nSolution found!\n");
    }
    else
    {
        printf("\nUnable to find a solution!\n");
    }
#ifdef WIN32
    int t2 = GetTickCount();
    printf("\nTime spent: %ims\n", t2-t1);
#endif
    return 0;
}
