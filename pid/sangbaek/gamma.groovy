package pid.sangbaek

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

class gamma{
  static def findGamma = { event ->
    def pbank = event.getBank("REC::Particle")
    return pbank -> (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==22}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
  }
}