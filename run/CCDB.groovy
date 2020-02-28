package run

import org.jlab.utils.groups.IndexedTable;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;


class CCDB {

    ccdb = new ConstantsManager();
    //ccdb.init(Arrays.asList(new String[]{"/daq/tt/fthodo","/calibration/eb/rf/config"}));
    ccdb.init(Arrays.asList(new String[]{"/daq/tt/fthodo","/calibration/eb/rf/config","/calibration/eb/rf/offset"}));
    rfTable = ccdb.getConstants(runNum,"/calibration/eb/rf/config");
    if (rfTable.hasEntry(1, 1, 1)){
            System.out.println(String.format("RF period from ccdb for run %d: %f",runNum,rfTable.getDoubleValue("clock",1,1,1)));
            rfPeriod = (float)rfTable.getDoubleValue("clock",1,1,1);
    }

}