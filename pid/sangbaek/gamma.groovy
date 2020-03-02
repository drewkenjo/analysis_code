package pid.sangbaek

import org.jlab.clas.physics.Vector3

class gamma{

  static def find_byEnergy = { gam_E ->
  	gam_E > 5
  }

  static def findGamma = { event ->
    def pbank = event.getBank("REC::Particle")
    return (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==22}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
  }
}