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

class DVCS {
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
      // if(twogamma.findsecondGamma(event) >0) return [null,null,null]
    }

    def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

    def parts = [11,2212,22].withIndex()
      .collect{pid,i -> new Particle(pid, *['px', 'py', 'pz'].collect{partbank.getFloat(it, inds[i])}) }

    def Vangle = {v1, v2 -> 
       if( v1.mag() * v2.mag() !=0 && v1.dot(v2)<v1.mag()*v2.mag() ) return Math.toDegrees( Math.acos(v1.dot(v2)/(v1.mag()*v2.mag()) ) ); 
    }

    // cut by kinematics
    // incoming
    def beam = new Particle(11, 0, 0, 10.6)// 10.6)
    def target = new Particle(2212, 0,0,0)

    // outgoing
    def (ele, pro, gam) = parts
    def W_vec = new Particle(beam)
    W_vec.combine(target, 1)
    W_vec.combine(ele, -1)
    def W = W_vec.mass()

    def GS = new Particle(beam)
    GS.combine(ele,-1)
    def VGS = GS.vector()
    def xB = -VGS.mass2()/(2*0.938*VGS.e())
    def Q2 = -VGS.mass2()
    if (W<2 || Q2<1 ||gam.e()<1) return [null, null, null]


    // cut by exclusivity
    def excl_cut = false

    def VG1 = gam.vector()
    def VE = ele.vector()
    def VPROT = pro.vector()

    def epX = new Particle(beam)
    epX.combine(target, 1)
    epX.combine(ele,-1)
    epX.combine(pro,-1)
    def VmissG = epX.vector()

    def egX = new Particle(beam)
    egX.combine(target, 1)
    egX.combine(ele,-1)
    egX.combine(gam,-1)
    def VmissP = egX.vector()

    def epgX = new Particle(beam)
    epgX.combine(target, 1)
    epgX.combine(ele,-1)
    epgX.combine(pro,-1)
    epgX.combine(gam,-1)
    def VMISS = epgX.vector()

    def Vhadr = (VPROT.vect()).cross(VGS.vect());
    def Vhad2 = (VGS.vect()).cross(VG1.vect());

    if( 
    // haz_g2==-1&& 
     VG1.e()>3
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
      excl_cut = true 
    }

    if (!excl_cut) return [null, null, null]


    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[inds[it]]]}
  }
}