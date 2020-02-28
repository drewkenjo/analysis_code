import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.sangbaek.DVCS

def hmm2_ep = new H1F("hmm2_ep", "missing mass squared, ep", 100,-2,4)
def hmm2_eg = new H1F("hmm2_eg", "missing mass squared, eg", 100,-2,4)
def hmm2_epg = new H1F("hmm2_epg", "missing mass squared, epg", 100,-0.2,0.2)
def hangle_epg = new H1F("hangle_epg", "Angle between gamma and epX", 100,-5 ,75)
def hangle_ep_eg = new H1F("hangle_ep_eg", "Angle between two planes, ep and eg", 190,-5,185)
def beam = new Particle(11, 0,0,10.6)//5)
def target = new Particle(2212, 0,0,0)

def h_kine_ele = new H2F("h_kine_ele", "e Kinematics", 100,0,40, 100, 0, 12)
def h_kine_pro = new H2F("h_kine_pro", "p Kinematics", 100,0,120, 100, 0, 12)
def h_kine_gam = new H2F("h_kine_gam", "#gamma Kinematics", 100,0,40, 100, 0, 12)
def h_Q2_xB = new H2F("h_Q2_xB", "Q^{2} - x_{B}",100,0,1,100,0,12);

def h_ele_rate = new H1F("h_ele_rate", "h_ele_rate",20,0,90)
def h_pro_rate = new H1F("h_pro_rate", "h_pro_rate",20,0,90)
def h_gam_rate = new H1F("h_gam_rate", "h_gam_rate",20,0,90)

def h_ep_azimuth = new H2F("h_ep_azimuth", "h_ep_azimuth",80,0,360,80, 0,360)
def h_ep_azimuth_diff = new H1F("h_ep_azimuth_diff", "h_ep_azimuth_diff",80,0,360)
def h_ep_polar = new H2F("h_ep_polar", "h_ep_polar",90,0,90,90,0,90)

def h_cross_section = new H1F("h_cross_section","DVCS cross section", 36,0,360)

def h_totalevent = new H1F("h_totalevent","total events",1,0,1)

def h_ele_phi = (1..6).collect{
  sec_num=it
  def hist = new H1F("electron phi sec"+sec_num,80,0,360)
  return hist
}
def h_W_sec = (1..6).collect{
  sec_num=it
  def hist = new H1F("h_W_"+sec_num,"h_W_"+sec_num,100,0,10)
  return hist
}

def h_Q2_xB_sec = (1..6).collect{
  sec_num=it
  def hist = new H2F("h_Q2_xB_"+sec_num,"h_Q2_xB_"+sec_num,100,0,1,100,0,12)
  return hist
}

def h_t_sec = (1..6).collect{
  sec_num=it
  def hist = new H1F("h_t_"+sec_num,"h_t_"+sec_num,100,0,4)
  return hist
}

def h_phi_sec = (1..6).collect{
  sec_num=it
  def hist = new H1F("h_phi"+sec_num,360,0,360)
  return hist
}

def h_y_sec = (1..6).collect{
  sec_num=it
  def hist = new H1F("W"+sec_num,100,0,1)
  return hist
}


for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)
h_totalevent.setBinContent(0,h_totalevent.getBinContent(0)+reader.getSize())


while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def dsets= DVCS.getEPG(event)
    def (ele, pro, gam) = dsets*.particle
    def (ele_sec, pro_sec, gam_sec) = dsets*.sector
    def Vangle = {v1, v2 -> 
       if( v1.mag() * v2.mag() !=0 && v1.dot(v2)<v1.mag()*v2.mag() ) return Math.toDegrees( Math.acos(v1.dot(v2)/(v1.mag()*v2.mag()) ) ); 
    }

    if(ele!=null) {
      def epX = new Particle(beam)
      epX.combine(target, 1)
      epX.combine(ele,-1)
      epX.combine(pro,-1)

      def egX = new Particle(beam)
      egX.combine(target, 1)
      egX.combine(ele,-1)
      egX.combine(gam,-1)

      def epgX = new Particle(beam)
      epgX.combine(target, 1)
      epgX.combine(ele,-1)
      epgX.combine(pro,-1)
      epgX.combine(gam,-1)

      def GS = new Particle(beam)
      GS.combine(ele,-1)

      def W = new Particle(GS)
      GS.combine(target,1)

      // def VG1 = gam.vector()
      def VGS = GS.vector()
      // def VMISS = epgX.vector()
      // def VmissP = egX.vector()
      // def VmissG = epX.vector()
      def VB = beam.vector()
      def VE = ele.vector()
      def Vlept = (VB.vect()).cross(VE.vect());

      def VPROT = pro.vector()
      def Vhadr = (VPROT.vect()).cross(VGS.vect());

      def TrentoAng = (float)Vangle(Vlept,Vhadr);

      def Q2 = -VGS.mass2()
      def xB = -VGS.mass2()/(2*0.938*VGS.e())

      def mand = new Particle(pro)
      mand.combine(target,-1)
      def Vmand = mand.vector()
      def t = -Vmand.mass2() //-t

      mom_gam = gam.vector().vect()
      mom_epX = epX.vector().vect()

      norm_ep = ele.vector().vect().cross(pro.vector().vect())
      norm_eg = ele.vector().vect().cross(gam.vector().vect())

      // if (Math.toDegrees(ele.phi())>0 && Math.toDegrees(ele.phi())<4.5){
      h_ele_rate.fill(Math.toDegrees(ele.theta()))
      // }
      // if (Math.toDegrees(pro.phi())>0 && Math.toDegrees(pro.phi())<4.5){
      h_pro_rate.fill(Math.toDegrees(pro.theta()))
      // }
      // if (Math.toDegrees(pro.phi())>0 && Math.toDegrees(gam.phi())<4.5){      
      h_gam_rate.fill(Math.toDegrees(gam.theta()))
      // }
      def ele_phi = Math.toDegrees(ele.phi())
      if (ele_phi<0) ele_phi=360+ele_phi
      def pro_phi = Math.toDegrees(pro.phi())
      if (pro_phi<0) pro_phi=360+pro_phi
      h_ep_azimuth.fill(pro_phi,ele_phi)
      h_ep_azimuth_diff.fill(Math.abs(pro_phi-ele_phi))
      h_ep_polar.fill(Math.toDegrees(pro.theta()),Math.toDegrees(ele.theta()))
      hmm2_ep.fill(epX.mass2())
      hmm2_eg.fill(egX.mass2())
      hmm2_epg.fill(epgX.mass2())
      hangle_epg.fill(Vangle(mom_gam,mom_epX))
      hangle_ep_eg.fill(Vangle(norm_ep,norm_eg))
      h_kine_ele.fill(Math.toDegrees(ele.vector().vect().theta()),ele.vector().vect().mag())
      h_kine_pro.fill(Math.toDegrees(pro.vector().vect().theta()),pro.vector().vect().mag())
      h_kine_gam.fill(Math.toDegrees(gam.vector().vect().theta()),gam.vector().vect().mag())
      h_Q2_xB.fill(xB,Q2);

      if((VPROT.vect()).dot(Vlept)<0)TrentoAng=-TrentoAng;
      if (TrentoAng<0) TrentoAng = 360+TrentoAng
      if (Q2>1 && Q2<5 && xB<0.5 && xB>0.2 && t<0.5 && t>0.2) h_cross_section.fill(TrentoAng)

      h_ele_phi[ele_sec-1].fill(ele_phi)

      h_Q2_xB_sec[ele_sec-1].fill(xB,Q2)
      h_W_sec[ele_sec-1].fill(W.mass())
    }
  }
}
reader.close()
}


// lumi = (double) 1.0558*0.0001
// xsec = (double) 9.0285*100000
// tot_rate = lumi * xsec
// phi_acceptance = (double) 4.5/360
// ratio = (double) 2000000/tot_rate
// ratio = (double) ratio/phi_acceptance
// h_ele_rate.divide(ratio)
// h_pro_rate.divide(ratio)
// h_gam_rate.divide(ratio)

def out = new TDirectory()
out.mkdir('/spec')
out.cd('/spec')
out.addDataSet(h_totalevent)

out.mkdir('/epg')
out.cd('/epg')
out.addDataSet(hmm2_ep)
out.addDataSet(hmm2_epg)
out.addDataSet(hmm2_eg)
out.addDataSet(hangle_epg)
out.addDataSet(hangle_ep_eg)
out.addDataSet(h_kine_ele)
out.addDataSet(h_kine_pro)
out.addDataSet(h_kine_gam)
out.addDataSet(h_Q2_xB)

out.mkdir('/rates')
out.cd('/rates')
out.addDataSet(h_ele_rate)
out.addDataSet(h_pro_rate)
out.addDataSet(h_gam_rate)

out.mkdir('/angular')
out.cd('/angular')
out.addDataSet(h_ep_azimuth)
out.addDataSet(h_ep_polar)
out.addDataSet(h_ep_azimuth_diff)
h_ele_phi.each{out.addDataSet(it)}

out.mkdir('/xsec')
out.cd('/xsec')
out.addDataSet(h_cross_section)

out.writeFile('dvcs_out.hipo')
