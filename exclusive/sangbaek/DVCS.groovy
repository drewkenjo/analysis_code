package exclusive.sangbaek

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.sangbaek.gamma
import pid.sangbaek.proton
import pid.sangbaek.electron
import org.jlab.clas.physics.LorentzVector;
import event.Event
import utils.KinTool

class DVCS {

  static def getEPG(event, electron_selector, proton_selector, gamma_selector) {

    def findElectron = { ev ->
      def electron_candidate = electron_selector.applyCuts_Brandon(ev)
      if(electron_candidate) electron_candidate
      .max{ind -> 
        if(event.px[ind]&& event.py[ind]&& event.pz[ind]) (new Vector3(*[event.px, event.py, event.pz].collect{it[ind]})).mag2()}
    }

    def findProton = { ev ->
      def proton_candidate = proton_selector.applyCuts_Stefan(ev)
      if(proton_candidate) proton_candidate
      .max{ind -> 
        if(event.px[ind]&& event.py[ind]&& event.pz[ind]) (new Vector3(*[event.px, event.py, event.pz].collect{it[ind]})).mag2()}
    }

    def findGamma = { ev ->
      def gamma_candidate = gamma_selector.applyCuts_Stefan(ev)
      if(gamma_candidate) gamma_candidate
      .max{ind -> 
        if(event.px[ind]&& event.py[ind]&& event.pz[ind]) (new Vector3(*[event.px, event.py, event.pz].collect{it[ind]})).mag2()}
    }


    def inds = []
    for(def findPart in [findElectron, findProton, findGamma]) {
      def ind = findPart(event)
      inds.add(ind)
      
      // cut by pid
      if(ind == null) return [null, null, null]
    }

    def status = (0..2).collect{it -> event.status[inds[it]]}

    def pcal_sectors = event.pcal_sector
    def ei_sectors = event.ecal_inner_sector
    def eo_sectors = event.ecal_outer_sector
    def ftof_sectors = event.tof

    def secs = []

    (0..<3).each{
      def pcal_sector = pcal_sectors[inds[it]]
      def ei_sector = ei_sectors[inds[it]]
      def eo_sector = eo_sectors[inds[it]]
      if(pcal_sector) secs.add(pcal_sector)
      else if(ei_sector) secs.add(ei_sector)
      else if(eo_sector) secs.add(eo_sector)
      else if (it==1 && ftof_sectors[inds[it]]) secs.add(ftof_sectors[inds[it]].sector?.find{true})
      else secs.add(null)
    }

    def vz_ele = event.vz[inds[0]]
    def vz_pro = event.vz[inds[1]]
    def p_pro  = event.p[inds[1]]
    def vzdiff = Math.abs(vz_ele-vz_pro)
    if (vzdiff>2.5+2.5/p_pro) return [null, null, null]

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> 
        if(event.px[inds[i]]&& event.py[inds[i]]&& event.pz[inds[i]]) new Particle(pid, *[event.px, event.py, event.pz].collect{it[inds[i]]})
    }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[it], status:status[it]]}
  }

  static def getEPG_EB(event) {

    def findElectron = { ev -> (0..<ev.npart).find{ev.pid[it]==11 && ev.status[it]<0} }
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

    def status = (0..2).collect{it -> event.status[inds[it]]}

    def pcal_sectors = event.pcal_sector
    def ei_sectors = event.ecal_inner_sector
    def eo_sectors = event.ecal_outer_sector
    def ftof_sectors = event.tof

    def secs = []

    (0..<3).each{
      def pcal_sector = pcal_sectors[inds[it]]
      def ei_sector = ei_sectors[inds[it]]
      def eo_sector = eo_sectors[inds[it]]
      if(pcal_sector) secs.add(pcal_sector)
      else if(ei_sector) secs.add(ei_sector)
      else if(eo_sector) secs.add(eo_sector)
      else if (it==1 && ftof_sectors[inds[it]]) secs.add(ftof_sectors[inds[it]].sector?.find{true})
      else secs.add(null)
    }

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *[event.px, event.py, event.pz].collect{it[inds[i]]})
    }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[it], status:status[it]]}
  }

  static def getEPG_MC(event) {

    def findElectron = { ev -> (0..<ev.mc_npart).find{ev.mc_pid[it]==11}}
    def findProton = {ev -> (0..<ev.mc_npart).find{ev.mc_pid[it]==2212}}
    def findGamma = {ev -> (0..<ev.mc_npart).find{ev.mc_pid[it]==22}}

    def inds = []

    for(def findPart in [findElectron, findProton, findGamma]) {
      def ind = findPart(event)
      inds.add(ind)
      // cut by pid
      if(ind == null) return [null, null, null]
      // if(twogamma.findsecondGamma(event) >0) return [null,null,null]
    }

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *[event.mc_px, event.mc_py, event.mc_pz].collect{it[inds[i]]})
    }

    def secs = parts.collect{
      def phi = Math.toDegrees(it.theta())
      phi += 20;
      if (phi<0) phi+=360;
      int sec = (int) phi/60;
      return sec+1
    }

    def status = parts.collect{
      def theta = Math.toDegrees(it.theta())
      if (theta>2.5 && theta<5) return 1500 
      else if (theta<35) return 2500
      else return 4500
    }


    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[it], status:status[it]]}
  }



  static def KineCuts(xB, Q2, W, VE, VG1){
      xB<1 && W>2 && Q2>1 && VE.e()>2.1 && VG1.e() >3
  }

  static def ExclCuts(VG1, VE, VMISS, VmissP, VmissG, Vhadr, Vhad2){
    // haz_g2==-1&& 
     (KinTool.Vangle(VG1.vect(),VE.vect())>4
     && VMISS.e()<1.5 
     && VMISS.mass2() <0.2 && VMISS.mass2() >-0.2 
     && VmissP.mass2() < 3 && VmissP.mass2() > -0.25
     && VmissG.mass2() < 1 && VmissG.mass2() > -1
     && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 0.3
     && KinTool.Vangle(VG1.vect(),VmissG.vect()) < 3
     && KinTool.Vangle(Vhad2,Vhadr) < 25)
  }
}