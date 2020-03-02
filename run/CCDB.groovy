package run

import org.jlab.utils.groups.IndexedTable;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.utils.groups.IndexedTable;


class CCDB {

	int run_number

	ConstantsManager ccmanager // use standard ConstantsManger if wanted.

	//some useful cases for reading ccdb cases. other examples can be added
	IndexedTable rfTableConfig, rfTableOffset, targetTable
	Float rfPeriod
	Float rfOffset1, rfOffset2
	Float targetPosition, targetLength


	//default constrcutor
	def CCDB(){
	}

	//optional constructor with run number. automatically saves RF and target info
	def CCDB(int runNum, int verboseLevel=0){
		this.run_number = runNum
		this.ccmanager = new ConstantsManager()
        this.ccmanager.init(["/geometry/target","/calibration/eb/rf/config","/calibration/eb/rf/offset"]);
		this.readRF(verboseLevel)
		this.readTarget(verboseLevel)
	}

	//reads RF table
	def readRF(int verboseLevel=0){
		this.rfTableConfig = this.ccmanager.getConstants(this.run_number,"/calibration/eb/rf/config");
		if (this.rfTableConfig.hasEntry(1, 1, 1)){
			this.rfPeriod = (float) this.rfTableConfig.getDoubleValue("clock",1,1,1);
			if(verboseLevel>0) println(String.format("RF period from ccdb for run %d: %f",this.run_number,this.rfPeriod));
		}

		this.rfTableOffset = this.ccmanager.getConstants(this.run_number,"/calibration/eb/rf/offset");
		if (this.rfTableOffset.hasEntry(1, 1, 1)){
			this.rfOffset1 = (float)rfTableOffset.getDoubleValue("offset",1,1,1);
			this.rfOffset2 = (float)rfTableOffset.getDoubleValue("offset",1,1,2);
			if(verboseLevel>0) println(String.format("RF1 offset from ccdb for run %d: %f",run_number,rfOffset1));
			if(verboseLevel>0) println(String.format("RF2 offset from ccdb for run %d: %f",run_number,rfOffset2));
		}
	}

	//reads target location
	def readTarget(int verboseLevel=0){
		this.targetTable = this.ccmanager.getConstants(this.run_number,"/geometry/target");
		if (this.targetTable.hasEntry(0, 0, 0)){
			this.targetPosition = (float) this.targetTable.getDoubleValue("position",0,0,0);
			this.targetLength = (float) this.targetTable.getDoubleValue("length",0,0,0);
			if(verboseLevel>0) println(String.format("Target Position %d, %f",this.run_number,this.targetTable.getDoubleValue("position",0,0,0)));
			if(verboseLevel>0) println(String.format("Target Length %d, %f",this.run_number,this.targetTable.getDoubleValue("length",0,0,0)));
		}
	}
}