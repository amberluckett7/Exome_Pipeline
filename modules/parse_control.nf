process CONTROL_PARSER {
    tag "$controlFile"
    cpus 1
    input:
      path controlFile

    output:
      path "parsed_lines.txt", emit: reads
      path "versions.yml",      emit: versions

    script:
    """
    python <<EOF
import sys
from collections import defaultdict

families     = []  
trios        = [] 
maleList     = []
femaleList   = []
fastqRecords = []
sampleSet    = set() 
                 

with open("${controlFile}", "r") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        tokens = line.split()
        if line.startswith("FAMILY"):
            families.append(tokens[1:])  
        elif line.startswith("TRIO"):
            trios.append(tokens[1:])     
        elif line.startswith("MALE"):
            maleList.extend(tokens[1:])
        elif line.startswith("FEMALE"):
            femaleList.extend(tokens[1:])
        elif line.startswith("FASTQ"):
            record = {
                "sample"  : tokens[1],
                "platform": tokens[2],
                "flowcell": tokens[3],
                "fastq1"  : tokens[4],
                "fastq2"  : tokens[5]
            }
            fastqRecords.append(record)
            sampleSet.add(tokens[1])
        elif line.startswith("SPRING"):
            record = {
                "sample"  : tokens[1],
                "platform": tokens[2],
                "flowcell": tokens[3],
                "fastq1"  : tokens[4].removesuffix(".spring"),
                "fastq2"  : tokens[4].removesuffix(".spring")
            }
            fastqRecords.append(record)
            sampleSet.add(tokens[1])

# ----------------------------

familyDict = {}
for famLine in families:
    famName = "-".join(famLine)
    samples = famLine[0:]
    familyDict.setdefault(famName, set()).update(samples)

trioDict = {}
for triLine in trios:
    trioName = "-".join(triLine)
    samples  = triLine[0:]
    trioDict.setdefault(trioName, set()).update(samples)

sampleLaneCounts = defaultdict(int)
for rec in fastqRecords:
    sampleLaneCounts[rec["sample"]] += 1

familySizes = { fname: len(sampleSet) for fname, sampleSet in familyDict.items() }

def getSex(sample):
    isMale   = (sample in maleList)
    isFemale = (sample in femaleList)
    if isMale and isFemale:
        return "conflict"
    elif isMale:
        return "male"
    elif isFemale:
        return "female"
    else:
        return "NA"

def getFamily(sample):
    for fname, members in familyDict.items():
        if sample in members:
            return fname
    return "NA"

def getTrio(sample):
    # Return the first trio name in which this sample appears, else "NA"
    for tname, members in trioDict.items():
        if sample in members:
            return tname
    return "NA"

with open("parsed_lines.txt", "w") as outFile:
    for rec in fastqRecords:
        s       = rec["sample"]
        sex     = getSex(s)
        fam     = getFamily(s)
        trio    = getTrio(s)
        laneCount = sampleLaneCounts[s]
        famCount  = familySizes.get(fam, 0) if fam != "NA" else 0

        line = (
            f"sample={s};"
            f"platform={rec['platform']};"
            f"flowcell={rec['flowcell']};"
            f"sex={sex};"
            f"family={fam};"
            f"trio={trio};"
            f"fastq1={rec['fastq1']};"
            f"fastq2={rec['fastq2']};"
            f"sampleLaneCount={laneCount};"
            f"familySampleCount={famCount}"
        )
        outFile.write(line + "\\n")

with open("versions.yml", "w") as vf:
    vf.write(f\"\"\"${task.process}:
  control_parser: 1.0
\"\"\")
EOF

    """
}
