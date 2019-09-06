package exclusive.sangbaek

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.proton.Proton
import pid.electron.Electron
import pid.sangbaek.gamma
import pid.sangbaek.electron

class EPG {
  static def getEPG(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def findElectron = { event -> 
      Electron.findElectron(event)
    }
    def findElectron_pid = {event ->
      electron.find_byEVENT(event)
    }
    def findProton = { event ->
      Proton.findProton(event)
    }
    def findGamma = { event -> 
      gamma.findGamma.(event)
    }

    def inds = []
    for(def findPart in [findElectron, findProton, findGamma]) {
      def ind = findPart(event)
      inds.add(ind)
      if(ind==null)
        return [null, null, null]
    }

    def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *['px', 'py', 'pz'].collect{partbank.getFloat(it, inds[i])}) }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[inds[it]]]}
  }
}

