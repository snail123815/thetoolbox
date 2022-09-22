from pathlib import Path

def findStartLine(csvFile: Path) -> int:
    with csvFile.open('r') as file:
        for i, line in enumerate(file):
            if not line.strip()[0] == '#':
                if line != '':
                    return i
        return 0
