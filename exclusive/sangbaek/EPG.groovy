package exclusive.sangbaek

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.proton.Proton
import pid.electron.Electron
import pid.sangbaek.gamma
import pid.sangbaek.electron
import org.jlab.clas.physics.LorentzVector;

class EPG {
  static def getEPG(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def findElectron = { ev -> 
      electron.find_byEVENT(ev)
    }
    def findProton = { ev ->
      Proton.findProton(ev)
    }
    def findGamma = { ev -> 
      gamma.findGamma(ev)
    }

    def inds = []
    for(def findPart in [findElectron, findProton, findGamma]) {
      def ind = findPart(event)
      inds.add(ind)
      
      // cut by pid
      if(ind==null) return [null, null, null]

      // cut by kinematics
      def beam = new Particle(11, 0,0,10.6)//7.546)
      def target = new Particle(2212, 0,0,0)
      def W_vec = new LorentzVector(0,0,0,0)
      W_vec.add(VB)
      W_vec.add(VT)
      W_vec.sub(VE)
      def W = W_vec.mass()

      def GS = new Particle(beam)
      GS.combine(ele,-1)
      def VGS = GS.vector()
      def xB = -VGS.mass2()/(2*0.938*VGS.e())
      def Q2 = -VGS.mass2()
      if (W<2 || Q2<1) return [null, null, null]

      // cut by exclusivity
       if( haz_g2==-1
       && VG1.e()>3
       && Vangle(VG1.vect(),VE.vect())>4
       && VMISS.e()<1.5 
       && VMISS.mass2() <0.2 && VMISS.mass2() >-0.2 
       && VmissP.mass2() < 3 && VmissP.mass2() > -0.25
       && VmissG.mass2() < 1 && VmissG.mass2() > -1
       && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 0.3
       && Vangle(VG1.vect(),VmissG.vect()) < 3
       && Vangle(Vhad2,Vhadr) < 25
       //&& Vangle(Vhad2,Vhadr) < 90
      ){

    }

    def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *['px', 'py', 'pz'].collect{partbank.getFloat(it, inds[i])}) }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[inds[it]]]}
  }
}

