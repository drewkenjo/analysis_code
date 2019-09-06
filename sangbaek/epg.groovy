import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.sangbaek.EPG

def hmm2_ep = new H1F("hmm2_ep", "missing mass squared, ep", 200,-2,4)
def hmm2_eg = new H1F("hmm2_eg", "missing mass squared, eg", 200,-2,4)
def hmm2_epg = new H1F("hmm2_epg", "missing mass squared, epg", 200,-2,4)
def hangle_epg = new H1F("hangle_epg", "Angle between gamma and epX", 200,-5 ,75)
def hangle_ep_eg = new H1F("hangle_ep_eg", "Angle between two planes, ep and eg", 380,-5,185)
def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)

while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def (ele, pro, gam) = EPG.getEPG(event)*.particle
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

      mom_gam = gam.vector().vect()
      mom_epX = epX.vector().vect()

      norm_ep = ele.vector().vect().cross(pro.vector().vect())
      norm_eg = ele.vector().vect().cross(gam.vector().vect())

      hmm2_ep.fill(epX.mass2())
      hmm2_eg.fill(egX.mass2())
      hmm2_epg.fill(epgX.mass2())
      hangle_epg.fill(Vangle(mom_gam,mom_epX))
      hangle_ep_eg.fill(Vangle(norm_ep,norm_eg))
    }
  }
}

reader.close()
}

def out = new TDirectory()
out.mkdir('/epg')
out.cd('/epg')
out.addDataSet(hmm2_ep)
out.addDataSet(hmm2_epg)
out.addDataSet(hmm2_eg)
out.addDataSet(hangle_epg)
out.addDataSet(hangle_ep_eg)
out.writeFile('epg_out.hipo')
