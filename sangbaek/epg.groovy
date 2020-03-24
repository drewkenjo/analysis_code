import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.EPG
import exclusive.sangbaek.DVCS
import utils.KinTool
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import pid.electron.ElectronSelector
import pid.sangbaek.electron
import run.Run

// default run_number
def run_number = -1000 

// Define default beam and target
def beam = new Particle(11, 0, 0, 10.604)
def VB = beam.vector()
def target = new Particle(2212, 0, 0, 0)

// Define default run package
def run = new Run()
def torus_scale = "-1"
def field_setting = ["-1" : "inbending", "1" :"outbending"]

// Define event counter to provide processing status
def event_count = 0
def dvcs_count = 0
def file_count = 0

// define a electron selector
def electron_selector = new electron()

// define ordinal number
def ordinal={
    (it % 100 == 11 || it % 100 == 12 || it % 100 == 13) ? it + "th" : it + ["th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th"][it % 10];
}

// define histograms to be drawn.
// missing mass
def hmm2_ep = new H1F("hmm2_ep", "missing mass squared, ep", 100,-2,4)
def hmm2_eg = new H1F("hmm2_eg", "missing mass squared, eg", 100,-2,4)
def hmm2_epg = new H1F("hmm2_epg", "missing mass squared, epg", 100,-0.2,0.2)
def hmm2_epg_dvcs = new H1F("hmm2_epg_dvcs", "missing mass squared, epg (DVCS)", 100,-0.2,0.2)

// angle between planes
def hangle_epg = new H1F("hangle_epg", "Angle between gamma and epX", 100,-5 ,75)
def hangle_ep_eg = new H1F("hangle_ep_eg", "Angle between two planes, ep and eg", 190,-5,185)

// kinematic variables (momentum vs theta)
def h_ele_mom_theta = new H2F("h_ele_mom_theta", "e momentum vs theta", 100,0,40, 100, 0, 12)
def h_pro_mom_theta = new H2F("h_pro_mom_theta", "p momentum vs theta", 100,0,120, 100, 0, 12)
def h_gam_mom_theta = new H2F("h_gam_mom_theta", "#gamma momentum vs theta", 100,0,40, 100, 0, 12)

// electron phi (sanity check)
def h_ele_phi = (1..6).collect{
  def hist = new H1F("electron phi sec"+it,80,0,360)
  return hist
}

def h_Q2_xB = new H2F("h_Q2_xB", "Q^{2} - x_{B}",100,0,1,100,0,12);
def h_t_xB = new H2F("h_t_xB", "-t - x_{B}",100,0,1,100,0,2);
def h_Q2_t = new H2F("h_Q2_t", "Q^{2} - -t",100,0,2,100,0,12);

def h_Q2_xB_cond = [:].withDefault{new H2F("h_Q2_xB_$it", "Q2 vs xB $it", 100,0,1,100,0,12)}

def h_Q2_theta = new H2F("h_Q2_theta", "Q^{2} - theta",100,0,45,100,0,12);

def h_ele_rate = new H1F("h_ele_rate", "h_ele_rate",90,0,90)
def h_pro_rate = new H1F("h_pro_rate", "h_pro_rate",90,0,90)
def h_gam_rate = new H1F("h_gam_rate", "h_gam_rate",90,0,90)

def h_ep_azimuth = new H2F("h_ep_azimuth", "h_ep_azimuth",80,0,360,80, 0,360)
def h_ep_azimuth_diff = new H1F("h_ep_azimuth_diff", "h_ep_azimuth_diff",80,0,360)
def h_ep_polar = new H2F("h_ep_polar", "h_ep_polar",90,0,90,90,0,90)

def h_cross_section = [:].withDefault{new H1F("h_cross_section_$it","DVCS cross section_$it", 24,0,360)}

// count total events collected
def h_totalevents = new H1F("h_totalevents","total events",1,0,1)
def h_dvcsevents = new H1F("h_dvcsevents","DVCS events",1,0,1)

// sector dependent 
def h_W_sec = (1..6).collect{
  def hist = new H1F("h_W_"+it,"h_W_"+it,100,0,10)
  return hist
}

def h_Q2_xB_sec = (1..6).collect{
  def hist = new H2F("h_Q2_xB_"+it,"h_Q2_xB_"+it,100,0,1,100,0,12)
  return hist
}

def h_t_sec = (1..6).collect{
  def hist = new H1F("h_t_"+it,"h_t_"+it,100,0,4)
  return hist
}

def h_phi_sec = (1..6).collect{
  def hist = new H1F("h_phi"+it,360,0,360)
  return hist
}

def h_y_sec = (1..6).collect{
  def hist = new H1F("h_y"+it,100,0,1)
  return hist
}

def xB_array = [0.1, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35, 0.38, 0.42, 0.58]
def t_array = [0.09, 0.13, 0.18, 0.23, 0.3, 0.39, 0.52, 0.72, 1.1, 2]
def binnumber = {xB, theta, t ->
  if (xB <xB_array[0] || xB >xB_array[11] || t <t_array[0] || t >t_array[9]) return null
  int xBQbin = 2*xB_array.findIndexOf{ xB < it}-1
  if (xBQbin==-3) xBQbin = 21
  if (xBQbin>1 && Math.toDegrees(theta)<15) xBQbin--

  int tbin = t_array.findIndexOf{ t < it} -1
  if (tbin==-1) tbin = 8
  return 21*tbin + xBQbin
}

// main loop
for(fname in args) {

  // std print out file count
  file_count++
  println("Reading " + ordinal(file_count) + " file...")

  // run number from file name
  def name = fname.split('/')[-1]
  def m = name =~ /\d{4,6}/
  def run_number_from_file = m[0].toInteger()

  // for each run, change beam energy.
  if (run_number != run_number_from_file){
    run_number = run_number_from_file
    run.run_number = run_number
    run.ReadDB()
    beam.setProperty("pz", run.rcdb.rcdb_dict["beam_energy"]/1000.0)
    VB = beam.vector()

    // if polarity is different, renew electron cut parameters.
    if (String.format("%.0f",run.rcdb.rcdb_dict["torus_scale"]) != torus_scale){
      torus_scale = String.format("%.0f",run.rcdb.rcdb_dict["torus_scale"])
      electron_selector.electron_candidate.ebeam = run.rcdb.rcdb_dict["beam_energy"]/1000.0
      electron_selector.electron_candidate.setElectronCutParameters(field_setting[torus_scale])
    }
  }

  def reader = new HipoDataSource()
  reader.open(fname)

  while(reader.hasEvent()) {
    def dataevent = reader.getNextEvent()
    def event = EventConverter.convert(dataevent)

    // count events so that users can know the program is running
    event_count++
    if (event_count%500000 == 0){
      println("Processing "+0.1*event_count.intdiv(500000)+" M-th event...")
    }

    if (event.npart>0) {

      // get epg coincidence, no exclusive cut applied. electron cut from Brandon's package
      def dsets = DVCS.getEPG(event, electron_selector)
      def (ele, pro, gam) = dsets*.particle

      // process only if there's a epg set in coincidence
      if(ele!=null) {

        // get sector
        def (ele_sec, pro_sec, gam_sec) = dsets*.sector
        
        // get 4 momentum variables!
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

        def GS = new Particle(beam)
        GS.combine(ele,-1)
        def VGS = GS.vector()

        def W_vec = new Particle(beam)
        W_vec.combine(target, 1)
        W_vec.combine(ele, -1)
        def W = W_vec.mass()

        def VG1 = gam.vector()
        // def VB = beam.vector() // VB already defined above
        def VE = ele.vector()
        def Vlept = (VB.vect()).cross(VE.vect());

        def VPROT = pro.vector()
        def Vhadr = (VPROT.vect()).cross(VGS.vect());
        def Vhad2 = (VGS.vect()).cross(VG1.vect());


        // Now kinematics used to cross sections
        // xB
        def xB = -VGS.mass2()/(2*0.938*VGS.e())
        // Q2
        def Q2 = -VGS.mass2()
        // phi
        def TrentoAng = (float) KinTool.Vangle(Vlept,Vhadr);
        if((VPROT.vect()).dot(Vlept)<0)TrentoAng=-TrentoAng;
        if (TrentoAng<0) TrentoAng = 360+TrentoAng
        // -t
        def mand = new Particle(pro)
        mand.combine(target,-1)
        def Vmand = mand.vector()
        def t = -Vmand.mass2() //-t

        // Fill Histogram
        h_ele_rate.fill(Math.toDegrees(ele.theta()))
        h_pro_rate.fill(Math.toDegrees(pro.theta()))
        h_gam_rate.fill(Math.toDegrees(gam.theta()))

        h_ele_mom_theta.fill(Math.toDegrees(ele.theta()),ele.p())
        h_pro_mom_theta.fill(Math.toDegrees(pro.theta()),pro.p())
        h_gam_mom_theta.fill(Math.toDegrees(gam.theta()),gam.p())

        // Trento like angle from ep and eg plane
        norm_ep = ele.vector().vect().cross(pro.vector().vect())
        norm_eg = ele.vector().vect().cross(gam.vector().vect())
        hangle_ep_eg.fill(KinTool.Vangle(norm_ep,norm_eg))

        // reconstructed and detected gamma angle deviation
        hangle_epg.fill(KinTool.Vangle(VG1.vect(),VmissG.vect()))

        // check deviation from elastic
        def ele_phi = Math.toDegrees(ele.phi())
        if (ele_phi<0) ele_phi=360+ele_phi
        def pro_phi = Math.toDegrees(pro.phi())
        if (pro_phi<0) pro_phi=360+pro_phi
        h_ep_azimuth.fill(pro_phi,ele_phi)
        h_ep_azimuth_diff.fill(Math.abs(pro_phi-ele_phi))
        h_ep_polar.fill(Math.toDegrees(pro.theta()),Math.toDegrees(ele.theta()))

        // check missing mass
        hmm2_ep.fill(epX.mass2())
        hmm2_eg.fill(egX.mass2())
        hmm2_epg.fill(epgX.mass2())


        // kinematic range
        h_Q2_xB.fill(xB,Q2);
        h_Q2_t.fill(t,Q2)
        h_t_xB.fill(xB,t)

        // working on binning
        double theta_e = Math.toDegrees(ele.theta()) /5
        String theta_label = String.format("%.0f_theta_%.0f",5*theta_e.trunc(),5*theta_e.trunc()+5)
        h_Q2_xB_cond[theta_label].fill(xB,Q2)
        if (W>2) h_Q2_xB_cond["W>2"].fill(xB,Q2)
        else h_Q2_xB_cond["W<2"].fill(xB,Q2)
        h_Q2_theta.fill(Math.toDegrees(ele.theta()),Q2);
        
        if (event.status[dsets.pindex[1]]>=4000) h_Q2_xB_cond['pro_CD'].fill(xB,Q2)
        if (event.status[dsets.pindex[2]]<2000) h_Q2_xB_cond['gam_FT'].fill(xB,Q2)

        h_ele_phi[ele_sec-1].fill(ele_phi)
        h_Q2_xB_sec[ele_sec-1].fill(xB,Q2)
        h_W_sec[ele_sec-1].fill(W)
        h_t_sec[ele_sec-1].fill(t)
        h_phi_sec[ele_sec-1].fill(TrentoAng) 
        h_y_sec[ele_sec-1].fill(KinTool.calcY(VB, VE))

        // exclusive cuts
        if (DVCS.ExclCuts(VG1, VE, VMISS, VmissP, VmissG, Vhadr, Vhad2)){
          dvcs_count++
          // if (Q2>1 && Q2<5 && xB<0.5 && xB>0.2 && t<0.5 && t>0.2) h_cross_section.fill(TrentoAng)
          def bin_number = binnumber(xB, ele.theta(), t)
          h_cross_section['dvcs_'+bin_number].fill(TrentoAng)
          if (event.status[dsets.pindex[1]]>=4000){
            h_Q2_xB_cond['dvcs_pro_CD_'+bin_number].fill(xB,Q2)
            h_cross_section['pro_CD_'+bin_number].fill(TrentoAng)
          }
          if (event.status[dsets.pindex[2]]<2000){
            h_Q2_xB_cond['dvcs_gam_FT_'+bin_number].fill(xB,Q2)
            h_cross_section['gam_FT_'+bin_number].fill(TrentoAng)            
          }
          if (event.status[dsets.pindex[1]]>=4000 && event.status[dsets.pindex[2]]<2000){
            h_Q2_xB_cond['dvcs_pro_CD_gam_FT_'+bin_number].fill(xB,Q2)
            h_cross_section['pro_CD_gam_FT_'+bin_number].fill(TrentoAng)            
          }            
        } // exclusivity cuts ended
        //add here for analysis

      }// event with e, p, g
    }// event with particles
  }// file loop
  reader.close()
}// file closed

// event number histograms
h_totalevents.setBinContent(0, event_count)  
h_dvcsevents.setBinContent(0, dvcs_count)  

// // If luminosities are calculated...
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
out.addDataSet(h_totalevents)
out.addDataSet(h_dvcsevents)

out.mkdir('/epg')
out.cd('/epg')
out.addDataSet(hmm2_ep)
out.addDataSet(hmm2_epg)
out.addDataSet(hmm2_eg)
out.addDataSet(hmm2_epg_dvcs)
out.addDataSet(hangle_epg)
out.addDataSet(hangle_ep_eg)
out.addDataSet(h_ele_mom_theta)
out.addDataSet(h_pro_mom_theta)
out.addDataSet(h_gam_mom_theta)
out.addDataSet(h_Q2_xB)
out.addDataSet(h_Q2_theta)

out.addDataSet(h_t_xB)
out.addDataSet(h_Q2_t)
h_Q2_xB_cond.values().each{out.addDataSet(it)}

h_Q2_xB_sec.each{out.addDataSet(it)}
h_W_sec.each{out.addDataSet(it)}
h_t_sec.each{out.addDataSet(it)}
h_phi_sec.each{out.addDataSet(it)}
h_y_sec.each{out.addDataSet(it)}


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
h_cross_section.values().each{out.addDataSet(it)}

out.writeFile('dvcs_out.hipo')