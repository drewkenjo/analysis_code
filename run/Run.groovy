package run
import run.CCDB
import run.RCDB

class Run {

    int run_number
    CCDB ccdb
    RCDB rcdb
    def special_runs =[:]

    // default constructor
    def Run(){
        this.run_number = -1000
    }

    // constructor with run number
    def Run(int runNum, int verboseLevel=0){
        this.run_number = runNum
        this.ReadDB(verboseLevel)
        this.specialRuns()
    }

    // read both DB at once
    def ReadDB(int verboseLevel=0){
        if(this.run_number != -1000){
            this.ReadCCDB(verboseLevel)
            this.ReadRCDB(verboseLevel)            
        }
    }

    // read RCDB
    def ReadCCDB(int verboseLevel=0){
        this.ccdb = new CCDB(this.run_number, verboseLevel)
    }

    def ReadRCDB(int verboseLevel=0){
        this.rcdb = new RCDB(this.run_number, verboseLevel)
    }

    def specialRuns(){
        List runs_fall2018Inbending = [4887,4888,4900,4893,4895,5418,5419]
        List runs_fall2018Outbending = [5443,5444,5453]
        this.special_runs["fall2018Inbending"] = runs_fall2018Inbending
        this.special_runs["fall2018Outbending"] = runs_fall2018Outbending
    }
}