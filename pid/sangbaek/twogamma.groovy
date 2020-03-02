package pid.sangbaek

import org.jlab.clas.physics.Vector3
import pid.sangbaek.gamma

class twogamma{

  static def find_byEnergy = { gam_E ->
  	gam_E > 5
  }

  static def findsecondGamma = { event ->

  	def ind_gam = gamma.findGamma(event)

    def pbank = event.getBank("REC::Particle")
    return (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==22 && it != ind_gam}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
  }
}