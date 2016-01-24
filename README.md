# UnderrailMutagenPuzzleSolver
A very dumb, but multi-threaded solver for hard Mutagen Puzzle in Underrail RPG

## Binaries

http://zxstudio.org/projects/underrail/UnderrailMutagenPuzzleSolver_1.0.zip

## Usage:

-f FILE -- read mutagens from file

-i      -- input mutagens by hand

-t N    -- use N threads (4 by default)

If -f is specified, the file should have this format:

MutagenName1: AA BB -CC DD EE -FF

MutagenName2: AA BB -CC DD EE -FF

The first mutagen in file is considered to be the goal.

## Example:

Exitus-1: WU JJ RJ LX RU IB LM RA D2 LS CI I5 DL IQ OY

Ovid-1: LX CW WU -RJ

Echo-2: P9 CI OY LS OC DL RJ -CW -IQ -WUetc...

Hint: to make search go faster, filter out impossible mutagens by hand,
for example those containing genes that are not in the final compound
and do not have a negative gene in any other sequence.

## Details:

This solver simply tries all possible combinations of mutagens until it finds one that works. I tried to squeeze a bit more performance out of it by making it multi-threaded and optimizing a few places, but it's still slow and dumb. It's only saving grace is that it works :)
