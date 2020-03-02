# CLAS12 RUN Class
This class is designed to sort each run's properties like polarities, luminosity, and beam energy. 

# Usage

Simple usage is as follows

```groovy
import run.Run

def run = new Run(5038)
println(run.rcdb.rcdb_dict)
println(run.special_runs)
```