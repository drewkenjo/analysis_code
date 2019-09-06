package pid.sangbaek

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

class electron{

  //1) default pid cut provided by EB (pid==11)
	static def find_byPID = { pid ->
    pid==11
  }

  //status check for trigger electron
  static def find_bySTATUS = { status ->
    status<0
  }

  //2) charge check
  static def find_byCHARGE = {charge ->
    charge<0
  }

  //3) NPhe cut
  static def find_byNPhe = {nphe ->
    nphe>2
  }

  //4) EC Inner vs EC outer cut
  // static def find_by_EC = { 

  // }

  //5) Sampling fraction cut
  static def find_bySampl = {sampl ->
    sampl >0.18
  }

  //6) PCAL Fiducial
  // static def find_by_PCAL = { 

  // }

  //7) DC R1, R2, R3 cut
  // static def find_by_DC = { 

  // }

  //8) Vertex-z positon cut
  static def find_byVZ = {vz, sec ->
    def vnom= -3 // for 5038
    def vstd = 5
    // switch(sec){
    //   case 1: vnom=-2.202; vstd=4.058
    //   case 2: vnom=-2.124; vstd=4.173
    //   case 3: vnom=-2.172; vstd=4.313
    //   case 4: vnom=-2.159; vstd=4.324
    //   case 5: vnom=-2.316; vstd=4.29
    //   case 6: vnom=-2.258; vstd=4.29
    // }
    Math.abs(vz-vnom) < 2.5*vstd
  }

  // FTOF Hit Response
  static def find_byFTOF = {evs, ind ->
    evs.getInt("pindex").each{pindex  ->  if (ind==pindex) return true}
    return false
  }

  // momentum cut
  static def find_byMOM = {mom, theta ->
    mom > 1.5 && theta>17*(1-mom/7)
  }

  static def find_byBANK = {pbank ->
    return (0..pbank.rows())
      .find{ ind->
        def pid = pbank.getInt('pid',ind)
        def status = pbank.getShort('status',ind)
        def charge = pbank.getByte('charge',ind)
        this.find_byPID(pid) && this.find_bySTATUS(status) && this.find_byCHARGE(charge)
      }
  }

  static def find_byEVENT = { event ->

    def pbank = event.getBank("REC::Particle")
    // Default pid, status, charge cut
    if (this.find_byBank(pbank)==null) return null
      // kinematics cut
    def ind = this.find_byBank(pbank)
    def mom = ['x','y','z'].collect{pbank.getFloat("p"+it,ind)}.sum()
    def pz = pbank.getFloat("pz", ind);
    def theta = Math.toDegrees(Math.acos(pz/mom))
    def vz = pbank.getFloat("vz",ind)
    def evc = event.getBank("REC::Calorimeter")
    e_ecal_E = 0

    //sampling fraction
    evc.getInt("pindex").eachWithIndex{ pindex, ind_c ->
      if (pindex==ind){
        det = evc.getInt("layer", ind_c);
        e_ecal_E+=evc.getFloat("energy",ind_c)
      }
    }
    sampl_frac = e_ecal_E/mom

    //nphe
    nphe = 0
    def evh = event.getBank("REC::Cherenkov")
    evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
      if(evh.getInt("detector",ind_h)==15 && pindex==ind) {
        nphe = evh.getFloat("nphe",ind_h)
      }
    }

    if(this.find_byMOM(mom,theta) && this.find_byVZ && this.find_bySamp(samp_frac) && this.find_byNPhe)) return ind
    else return null
  }
}