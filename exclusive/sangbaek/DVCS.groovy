package exclusive.sangbaek

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.proton.Proton
import pid.electron.Electron
import pid.sangbaek.gamma
import pid.sangbaek.twogamma
import pid.sangbaek.electron
import org.jlab.clas.physics.LorentzVector;
import event.Event
import utils.KinTool

class DVCS {

  static def getEPG(event, electron_ind) {

    def findElectron = { ev ->
      def electron_candidate = electron_ind.applyCuts_Brandon(ev)
      if(electron_candidate) electron_candidate
      .max{ind -> (new Vector3(*[event.px, event.py, event.pz].collect{it[ind]})).mag2()}
    }
    def findProton = {ev -> (0..<ev.npart).findAll{ev.pid[it]==2212}
      .max{ind -> (new Vector3(*[ev.px, ev.py, ev.pz].collect{it[ind]})).mag2()}
    }

    def findGamma = {ev -> (0..<ev.npart).findAll{ev.pid[it]==22}
      .max{ind -> (new Vector3(*[ev.px, ev.py, ev.pz].collect{it[ind]})).mag2()}
    }


    def inds = []
    for(def findPart in [findElectron, findProton, findGamma]) {
      def ind = findPart(event)
      inds.add(ind)
      
      // cut by pid
      if(ind == null) return [null, null, null]
      // if(twogamma.findsecondGamma(event) >0) return [null,null,null]
    }

    def secs1 = event.pcal_sector
    def secs2 = event.ecal_inner_sector
    def secs3 = event.ecal_outer_sector

    def secs = []

    (0..<3).each{
      def sec1 = secs1[inds[it]]
      def sec2 = secs2[inds[it]]
      def sec3 = secs3[inds[it]]
      if(sec1) secs.add(sec1)
      else if(sec2) secs.add(sec2)
      else if(sec3) secs.add(sec3)
      else secs.add(null)
    }

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *[event.px, event.py, event.pz].collect{it[inds[i]]})
    }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[it]]}
  }

  static def KineCuts(Q2, W, VG1){
      W>2 && Q2>1 && VG1.e() >1
  }

  static def ExclCuts(VG1, VE, VMISS, VmissP, VmissG, Vhadr, Vhad2){
    // haz_g2==-1&& 
     (VG1.e() > 3 
     && KinTool.Vangle(VG1.vect(),VE.vect())>4
     && VMISS.e()<1.5 
     && VMISS.mass2() <0.2 && VMISS.mass2() >-0.2 
     && VmissP.mass2() < 3 && VmissP.mass2() > -0.25
     && VmissG.mass2() < 1 && VmissG.mass2() > -1
     && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 0.3
     && KinTool.Vangle(VG1.vect(),VmissG.vect()) < 3
     && KinTool.Vangle(Vhad2,Vhadr) < 25
     && KinTool.Vangle(Vhad2,Vhadr) < 90)
  }
}