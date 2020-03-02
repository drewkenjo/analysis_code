package run

import org.jlab.detector.calib.utils.RCDBProvider
import org.jlab.detector.calib.utils.RCDBConstants

class RCDB {

    //https://github.com/JeffersonLab/clas12-offline-software/tree/development/common-tools/clas-detector/src/main/java/org/jlab/detector/calib/utils
    RCDBProvider provider // use standard rcdb provider if wanted.
    RCDBConstants constants // use standard rcdb constants if wanted.

    //run number to read rcdb
    int run_number

    // Keys for rcdb
    def rcdb_keys=[]
    def rcdb_dict =[:]
    // default constructor
    def RCDB(){
        this.run_number = -1000
        this.defineKeys()
    }

    // optional constructor with run number
    def RCDB(int runNum, int verboseLevel=0){
        this.run_number = runNum
        this.defineKeys()
        this.readRCDB(verboseLevel)
    }

    def defineKeys(){
        this.rcdb_keys << "beam_current"
        this.rcdb_keys << "beam_current_request"
        this.rcdb_keys << "beam_energy"
        this.rcdb_keys << "event_count"
        this.rcdb_keys << "events_rate"
        this.rcdb_keys << "evio_files_count"
        this.rcdb_keys << "half_wave_plate"
        // this.rcdb_keys << "is_valid_run_end"
        // this.rcdb_keys << "json_cnd"   
        this.rcdb_keys << "megabyte_count"
        this.rcdb_keys << "operators"
        this.rcdb_keys << "run_config"
        this.rcdb_keys << "run_end_time"
        this.rcdb_keys << "run_start_time"
        this.rcdb_keys << "run_type"
        this.rcdb_keys << "solenoid_current"
        this.rcdb_keys << "solenoid_scale"
        this.rcdb_keys << "status"
        this.rcdb_keys << "target"
        this.rcdb_keys << "torus_current"
        this.rcdb_keys << "torus_scale"
        this.rcdb_keys << "user_comment"
    }

    def readRCDB(int verboseLevel=0){
        if(this.run_number == -1000){
            println("Please provide run number. Usage:")
            println("def rcdb = new RCDB(5038)\n, or,")
            println("def rcdb = new RCDB()\n int run = 5038 \n rcdb.readRCDB()")
            return
        }
        else{
            this.provider = new RCDBProvider()
            this.constants = this.provider.getConstants(this.run_number)
            if(verboseLevel>0)  this.constants.show()
            this.rcdb_keys.each{key->
                // if(this.constants.getDouble(key)!=null){
                //     this.rcdb[key] = this.constants.getDouble(key)
                // }
                // if(this.constants.getString(key)!=null){
                //     this.rcdb[key] = this.constants.getString(key)
                // }
                // if(this.constants.getLong(key)!=null){
                //     this.rcdb[key] = this.constants.getLong(key)
                // }
                // if((Time)this.constants.getTime(key)!=null){
                //     this.rcdb[key] = this.constants.getTime(key)
                // }
                if(this.constants.get(key)!=null){
                    this.rcdb_dict[key] = this.constants.get(key).getValue()
                }
                else this.rcdb_dict[key] = null
            }
        }
    }
}