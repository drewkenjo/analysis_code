package sangbaek
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
import java.util.concurrent.ConcurrentHashMap


class dvcs{

  //defining histograms
  def hists = new ConcurrentHashMap()

  // missing mass
  def hmm2 = {new H1F("$it", "$it", 100, -2, 4)}

  // angle between planes
  def h_angle = {new H1F("$it", "$it", 190, -5 ,185)}

  // kinematic variables (momentum vs theta)
  def h_mom_theta = {new H2F("$it", "$it", 300, 0, 120, 100, 0, 12)}

  def h_Q2_xB = {new H2F("$it", "$it", 100, 0, 1,100, 0, 12)}
  def h_t_xB = {new H2F("$it", "$it", 100, 0, 1,100, 0, 2)}
  def h_Q2_t = {new H2F("$it", "$it", 100, 0, 2,100, 0, 12)}
  def h_Q2_theta = {new H2F("$it", "$it", 100, 0, 45, 100, 0, 12)}

  // polar angle
  def h_polar_rate = {new H1F("$it", "$it", 90, 0, 90)}
  // electron phi (sanity check)
  def h_azimuth_rate = {new H1F("$it","$it", 80, 0, 360)}


  def h_ep_azimuth = {new H2F("$it", "$it", 80, 0, 360, 80, 0, 360)}
  def h_ep_azimuth_diff = {new H1F("$it", "$it", 80, 0, 360)}
  def h_ep_polar = {new H2F("$it", "$it", 90, 0, 90, 90, 0, 90)}

  def h_cross_section = {new H1F("$it","$it", 24, 0, 360)}

  // count total events collected
  def h_events = {new H1F("$it","$it",10, 0,10)}

  // sector dependence of kinematic variables 
  def h_W = {new H1F("$it","$it",100,0,10)}
  def h_t = {new H1F("$it","$it",100,0,4)}
  def h_y = {new H1F("$it",100,0,1)}

  //binning
  def xB_array = [0.1, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35, 0.38, 0.42, 0.58]
  def t_array = [0.09, 0.13, 0.18, 0.23, 0.3, 0.39, 0.52, 0.72, 1.1, 2]
  def binnumber = {xB, theta, t ->
    if (xB <xB_array[0] || xB >xB_array[11] || t <t_array[0] || t >t_array[9]) return null
    int xBQbin = 2*xB_array.findIndexOf{ xB < it}-1
    if (xBQbin==-3) xBQbin = 21
    if (xBQbin>1 && Math.toDegrees(theta)<10) xBQbin--

    int tbin = t_array.findIndexOf{ t < it} -1
    if (tbin==-1) tbin = 8
    return 21*tbin + xBQbin
  }

  def electron_selector = new electron()
  def beam = LorentzVector.withPID(11, 0, 0, 10.6)
  def target = LorentzVector.withPID(2212, 0, 0, 0)

  def processEvent(event){

    if (event.npart>0) {
      
      // get epg coincidence, no exclusive cut applied. electron cut from Brandon's package
      def dsets = DVCS.getEPG(event, electron_selector)
      def (ele, pro, gam) = dsets*.particle.collect{it ? it.vector() : null} 
      // process only if there's a epg set in coincidence
      if(ele!=null) {
        // event number histograms
        hists.computeIfAbsent("/events/events", h_events).fill(0.5)  

        // get sector
        def (ele_sec, pro_sec, gam_sec) = dsets*.sector
        
        // get 4 momentum variables!
        def VmissG = beam + target - ele - pro
        def VmissP = beam + target - ele - gam
        def VMISS = beam + target - ele - pro -gam
        def VGS = beam - ele 
        def W = (beam + target -ele).mass()
        def Vlept = (beam.vect()).cross(ele.vect());
        def Vhadr = (pro.vect()).cross(VGS.vect());
        def Vhad2 = (VGS.vect()).cross(gam.vect());


        // Now kinematics used to cross sections
        def xB = KinTool.calcXb(beam, ele)
        def Q2 = KinTool.calcQ2(beam, ele)
        def TrentoAng = KinTool.calcPhiTrento(beam, ele, pro); // phi
        def t = KinTool.calcT(pro) //-t

        // Fill Histogram
        hists.computeIfAbsent("/epg/elec_polar_sec"+ele_sec, h_polar_rate).fill(Math.toDegrees(ele.theta()))
        hists.computeIfAbsent("/epg/prot_polar", h_polar_rate).fill(Math.toDegrees(pro.theta()))
        hists.computeIfAbsent("/epg/gam_polar", h_polar_rate).fill(Math.toDegrees(gam.theta()))

        hists.computeIfAbsent("/epg/elec_mom_theta_sec"+ele_sec, h_mom_theta).fill(Math.toDegrees(ele.theta()),ele.p())
        hists.computeIfAbsent("/epg/prot_mom_theta", h_mom_theta).fill(Math.toDegrees(pro.theta()),pro.p())
        hists.computeIfAbsent("/epg/gam_mom_theta", h_mom_theta).fill(Math.toDegrees(gam.theta()),gam.p())

        // Trento like angle from ep and eg plane
        hists.computeIfAbsent("/epg/hangle_ep_ep", h_angle).fill(KinTool.Vangle(ele.vect().cross(pro.vect()), ele.vect().cross(gam.vect())))
        // reconstructed and detected gamma angle deviation
        hists.computeIfAbsent("/epg/hangle_epg", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))

        // check deviation from elastic
        def ele_phi = Math.toDegrees(ele.phi())
        if (ele_phi<0) ele_phi=360+ele_phi
        def pro_phi = Math.toDegrees(pro.phi())
        if (pro_phi<0) pro_phi=360+pro_phi
        def gam_phi = Math.toDegrees(gam.phi())
        if (gam_phi<0) gam_phi=360+gam_phi

        hists.computeIfAbsent("/epg/elec_azimuth_sec"+ele_sec, h_azimuth_rate).fill(ele_phi)
        hists.computeIfAbsent("/epg/prot_azimuth", h_azimuth_rate).fill(pro_phi)
        hists.computeIfAbsent("/epg/gam_azimuth", h_azimuth_rate).fill(gam_phi)
        
        hists.computeIfAbsent("/epg/h_ep_azimuth", h_ep_azimuth).fill(pro_phi,ele_phi)
        hists.computeIfAbsent("/epg/h_ep_azimuth_diff", h_ep_azimuth_diff).fill(Math.abs(pro_phi-ele_phi))
        hists.computeIfAbsent("/epg/h_ep_polar", h_ep_polar).fill(Math.toDegrees(pro.theta()),Math.toDegrees(ele.theta()))

        // check missing mass
        hists.computeIfAbsent("/epg/hmm2_ep", hmm2).fill(VmissG.mass2())
        hists.computeIfAbsent("/epg/hmm2_eg", hmm2).fill(VmissP.mass2())
        hists.computeIfAbsent("/epg/hmm2_epg", hmm2).fill(VMISS.mass2())


        // kinematic range
        hists.computeIfAbsent("/epg/h_Q2_xB", h_Q2_xB).fill(xB,Q2)
        hists.computeIfAbsent("/epg/h_Q2_t", h_Q2_t).fill(t,Q2)
        hists.computeIfAbsent("/epg/h_t_xB", h_t_xB).fill(xB,t)
        hists.computeIfAbsent("/epg/h_Q2_theta", h_Q2_theta).fill(Math.toDegrees(ele.theta()),Q2);

        // working on binning // don't have to fill this because binning is over.
        // double theta_e = Math.toDegrees(ele.theta()) /5
        // String theta_label = String.format("%.0f_theta_%.0f",5*theta_e.trunc(),5*theta_e.trunc()+5)
        // h_Q2_xB_cond[theta_label].fill(xB,Q2)
        if (W>2)  hists.computeIfAbsent("/epg/h_Q2_xB_W>2", h_Q2_xB).fill(xB,Q2)
        else hists.computeIfAbsent("/epg/h_Q2_xB_W<2", h_Q2_xB).fill(xB,Q2)
        

        if (event.status[dsets.pindex[1]]>=4000) hists.computeIfAbsent("/epg/h_t_xB_pro_CD", h_Q2_xB).fill(xB,Q2)
        if (event.status[dsets.pindex[2]]<2000) hists.computeIfAbsent("/epg/h_t_xB_gam_FT", h_Q2_xB).fill(xB,Q2)


        // check CD alignment
        if (event.status[dsets.pindex[1]]>=4000){
          hists.computeIfAbsent("/epg/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/epg/prot_azimuth_CD", h_azimuth_rate).fill(pro_phi)
          if (W>2){
            hists.computeIfAbsent("/epg/prot_polar_CD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/epg/prot_azimuth_CD_W>2", h_azimuth_rate).fill(pro_phi)
          }
        }
        else if (event.status[dsets.pindex[1]]<4000){
          hists.computeIfAbsent("/epg/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/epg/prot_azimuth_FD", h_azimuth_rate).fill(pro_phi)
          if (W>2){
            hists.computeIfAbsent("/epg/prot_polar_FD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/epg/prot_azimuth_FD_W>2", h_azimuth_rate).fill(pro_phi)
          }
        }



        hists.computeIfAbsent("/epg/h_Q2_xB_sec"+ele_sec, h_Q2_xB).fill(xB,Q2)
        hists.computeIfAbsent("/epg/h_W_sec"+ele_sec, h_W).fill(W)
        hists.computeIfAbsent("/epg/h_t_sec"+ele_sec, h_t).fill(t)
        hists.computeIfAbsent("/epg/h_phi_sec"+ele_sec, h_cross_section).fill(TrentoAng) 
        hists.computeIfAbsent("/epg/h_y_sec"+ele_sec, h_y).fill(KinTool.calcY(beam, ele))

        // exclusive cuts
        if (DVCS.KineCuts(Q2, W, gam) && DVCS.ExclCuts(gam, ele, VMISS, VmissP, VmissG, Vhadr, Vhad2)){
          // if (Q2>1 && Q2<5 && xB<0.5 && xB>0.2 && t<0.5 && t>0.2) h_cross_section.fill(TrentoAng)
          def bin_number = binnumber(xB, ele.theta(), t)
          hists.computeIfAbsent("/dvcs/h_phi_bin_$bin_number", h_cross_section).fill(TrentoAng)
          hists.computeIfAbsent("/events/events", h_events).fill(1.5)  
          
          if (event.status[dsets.pindex[1]]>=4000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_pro_CD_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_pro_CD_bin_$bin_number", h_cross_section).fill(TrentoAng)
            hists.computeIfAbsent("/events/events", h_events).fill(2.5)  
          }

          if (event.status[dsets.pindex[2]]<2000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_gam_FT_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_gam_FT_bin_$bin_number", h_cross_section).fill(TrentoAng)
            hists.computeIfAbsent("/events/events", h_events).fill(3.5)  
          }

          if (event.status[dsets.pindex[1]]>=4000 && event.status[dsets.pindex[2]]<2000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_pro_CD_gam_FT_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_pro_CD_gam_FT_bin_$bin_number", h_cross_section).fill(TrentoAng)
            hists.computeIfAbsent("/dvcs/hmm2_epg", hmm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/events/events", h_events).fill(4.5)  
          }            

          if (event.status[dsets.pindex[1]]>=4000){
            hists.computeIfAbsent("/dvcs/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/dvcs/prot_azimuth_CD", h_azimuth_rate).fill(pro_phi)
            if (W>2){
                hists.computeIfAbsent("/dvcs/prot_polar_CD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
                hists.computeIfAbsent("/dvcs/prot_azimuth_CD_W>2", h_azimuth_rate).fill(pro_phi)
            }
          }
          else if (event.status[dsets.pindex[1]]<4000){
            hists.computeIfAbsent("/dvcs/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/dvcs/prot_azimuth_FD", h_azimuth_rate).fill(pro_phi)
            if (W>2){
              hists.computeIfAbsent("/dvcs/prot_polar_FD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
              hists.computeIfAbsent("/dvcs/prot_azimuth_FD_W>2", h_azimuth_rate).fill(pro_phi)
            }
          }
        } // exclusivity cuts ended
        //add here for analysis
      }// event with e, p, g
    }// event with particles
  }
}