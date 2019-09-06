package exclusive.sangbaek

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.sangbaek.electron

class EPG {
  static def getEPG(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def findElectron = { event -> 
      def pbank = event.getBank("REC::Particle")
      return(0..<pbank.rows()).find{pbank.getInt('pid',it)==11 && pbank.getShort('status',it)<0}
    }
    def findElectron_pid = {event -> electron.find_byEVENT(event)}
    def findProton = { event ->
      def pbank = event.getBank(("REC::Particle")
      return (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==2212}
        .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
    }
    def findGamma = { event -> 
      def pbank = event.getBank(("REC::Particle")
      (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==22}
        .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
    }

    def inds = []
    for(def findPart in [findElectron_pid, findProton, findGamma]) {
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

