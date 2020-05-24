package sangbaek

import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import exclusive.sangbaek.DVCS
import utils.KinTool
import event.Event
import event.EventConverter
import pid.electron.ElectronSelector
import pid.sangbaek.electron
import pid.sangbaek.proton
import pid.sangbaek.gamma
import run.Run
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase


class dvcs{

  //defining histograms
  def hists = new ConcurrentHashMap()

  // missing mass
  def hmm2 = {new H1F("$it", "$it", 100, -2, 4)}
  def hmm2_2 = {new H1F("$it", "$it", 100, -0.2, 0.2)}

  // angle between planes
  def h_angle = {new H1F("$it", "$it", 190, -5 ,185)}

  // kinematic variables (correlation)
  def h_theta_mom = {new H2F("$it", "$it", 120, 0, 12, 100, 0, 100)}
  def h_phi_mom = {new H2F("$it", "$it", 120, 0, 12, 360, -180, 180)}
  def h_theta_phi = {new H2F("$it", "$it", 360, -180, 180, 100, 0, 100)}
  def h_theta_t = {new H2F("$it", "$it", 100, 0, 4, 100, 0, 100)}
  def h_phi_t = {new H2F("$it", "$it", 100, 0, 4, 360, -180, 180)}
  def h_theta_trento = {new H2F("$it", "$it", 360, 0, 360, 100, 0, 100)}
  def h_phi_trento = {new H2F("$it", "$it", 360, 0, 360, 360, -180, 180)}

  def h_Q2_xB = {new H2F("$it", "$it", 100, 0, 1,100, 0, 12)}
  def h_t_xB = {new H2F("$it", "$it", 100, 0, 1,100, 0, 2)}
  def h_Q2_t = {new H2F("$it", "$it", 100, 0, 4,100, 0, 12)}
  def h_Q2_theta = {new H2F("$it", "$it", 100, 0, 45, 100, 0, 12)}

  // polar angle
  def h_polar_rate = {new H1F("$it", "$it", 360, 0, 90)}
  // electron phi (sanity check)
  def h_azimuth_rate = {new H1F("$it","$it", 80, 0, 360)}


  def h_ep_azimuth = {new H2F("$it", "$it", 80, 0, 360, 80, 0, 360)}
  def h_ep_azimuth_diff = {new H1F("$it", "$it", 80, 0, 360)}
  def h_ep_polar = {new H2F("$it", "$it", 90, 0, 90, 90, 0, 90)}

  def h_cross_section = {new H1F("$it","$it", 24, 0, 360)}

  // count total events collected
  def h_events = {new H1F("$it","$it",12, 0,12)}

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
    if (xBQbin>1 && Math.toDegrees(theta)<15) xBQbin--

    int tbin = t_array.findIndexOf{ t < it} -1
    if (tbin==-1) tbin = 8
    return 21*tbin + xBQbin
  }

  def phi_convention = {phi ->
    phi=phi + 180
    if (phi>180) phi = phi - 360
    if (phi<-180) phi = phi + 360
    return phi
  }

  def electron_selector = new electron()
  def proton_selector = new proton()
  def gamma_selector = new gamma()
  def beam = LorentzVector.withPID(11, 0, 0, 10.6)
  def target = LorentzVector.withPID(2212, 0, 0, 0)

  def processEvent(event){

    hists.computeIfAbsent("/events/events", h_events).fill(0.5)  

    if (event.npart>0) {

      hists.computeIfAbsent("/events/events", h_events).fill(1.5)  

      (0..<event.npart).findAll{event.pid[it]==2212}.each{ind->
        def prot = new Vector3(*[event.px, event.py, event.pz].collect{it[ind]})
        def prot_phi = Math.toDegrees(prot.phi())
        if (prot_phi<0) prot_phi=360+prot_phi
        if (event.status[ind]>=4000){
          hists.computeIfAbsent("/prot/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(prot.theta()))
          hists.computeIfAbsent("/prot/prot_azimuth_CD", h_azimuth_rate).fill(prot_phi)
        }
        else if (event.status[ind]<4000){
          hists.computeIfAbsent("/prot/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(prot.theta()))
          hists.computeIfAbsent("/prot/prot_azimuth_FD", h_azimuth_rate).fill(prot_phi)
        }
      }      
      // get epg coincidence, no exclusive cut applied. electron cut from Brandon's package
      def dsets = DVCS.getEPG(event, electron_selector, proton_selector, gamma_selector)
      def (ele, pro, gam) = dsets*.particle.collect{it ? it.vector() : null} 
      // process only if there's a epg set in coincidence
      if(ele!=null) {
        // event number histograms
        hists.computeIfAbsent("/events/events", h_events).fill(2.5)  

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

        hists.computeIfAbsent("/epg/elec_theta_mom_sec"+ele_sec, h_theta_mom).fill(ele.p(), Math.toDegrees(ele.theta()))
        hists.computeIfAbsent("/epg/prot_theta_mom", h_theta_mom).fill(pro.p(), Math.toDegrees(pro.theta()))
        hists.computeIfAbsent("/epg/gam_theta_mom", h_theta_mom).fill(gam.p(), Math.toDegrees(gam.theta()))

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
        hists.computeIfAbsent("/epg/hmm2_epg", hmm2_2).fill(VMISS.mass2())


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

        //calc tcol tmin
        def E = 10.6
        def M = PDGDatabase.getParticleMass(2212)
        def tmin = M*M*xB*xB/(1-xB+xB*M*M/Q2)
        def tcol = Q2*(Q2-2*xB*M*E)/xB/(Q2-2*M*E)
        // fill t dependence on 2 fold binning (xB, Q2)
        int xBbin = 1 + 2 * Math.floor(xB/0.2)
        int Q2bin = 1 + 2 * Math.floor(Q2/2)

        // hists.computeIfAbsent("/epg/corr/tmin",h_Q2_xB).fill(xB,Q2,tmin)
        // hists.computeIfAbsent("/epg/corr/tcol",h_Q2_xB).fill(xB,Q2,tcol)
        // hists.computeIfAbsent("/epg/corr/prot_theta_mom_xB_${xBbin}_Q2_${Q2bin}", h_theta_mom).fill(pro.p(), Math.toDegrees(pro.theta()))
        // hists.computeIfAbsent("/epg/corr/prot_phi_mom_xB_${xBbin}_Q2_${Q2bin}", h_phi_mom).fill(pro.p(), Math.toDegrees(pro.phi()-ele.phi()))
        // hists.computeIfAbsent("/epg/corr/prot_theta_phi_xB_${xBbin}_Q2_${Q2bin}", h_theta_phi).fill(Math.toDegrees(pro.phi()-ele.phi()), Math.toDegrees(pro.theta()))
        // hists.computeIfAbsent("/epg/corr/gam_phi_mom_xB_${xBbin}_Q2_${Q2bin}", h_phi_mom).fill(gam.p(), Math.toDegrees(gam.phi()-ele.phi()))
        // hists.computeIfAbsent("/epg/corr/gam_theta_mom_xB_${xBbin}_Q2_${Q2bin}", h_theta_mom).fill(gam.p(), Math.toDegrees(gam.theta()))
        // hists.computeIfAbsent("/epg/corr/gam_theta_phi_xB_${xBbin}_Q2_${Q2bin}", h_theta_phi).fill(Math.toDegrees(gam.phi()-ele.phi()), Math.toDegrees(gam.theta()))

        // exclusive cuts
        if (DVCS.KineCuts(xB, Q2, W, gam) && DVCS.ExclCuts(gam, ele, VMISS, VmissP, VmissG, Vhadr, Vhad2)){

          hists.computeIfAbsent("/dvcs/corr/tmin", h_Q2_xB).fill(xB,Q2,tmin)
          hists.computeIfAbsent("/dvcs/corr/tcol", h_Q2_xB).fill(xB,Q2,tcol)

          def pro_phi_convention = phi_convention(Math.toDegrees(pro.phi()-ele.phi()))
          def gam_phi_convention = phi_convention(Math.toDegrees(pro.phi()-ele.phi()))

          hists.computeIfAbsent("/dvcs/corr/prot_theta_mom_xB_${xBbin}_Q2_${Q2bin}", h_theta_mom).fill(pro.p(), Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/dvcs/corr/prot_phi_mom_xB_${xBbin}_Q2_${Q2bin}", h_phi_mom).fill(pro.p(), pro_phi_convention)
          hists.computeIfAbsent("/dvcs/corr/prot_theta_phi_xB_${xBbin}_Q2_${Q2bin}", h_theta_phi).fill(pro_phi_convention, Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/dvcs/corr/gam_phi_mom_xB_${xBbin}_Q2_${Q2bin}", h_phi_mom).fill(gam.p(), gam_phi_convention)
          hists.computeIfAbsent("/dvcs/corr/gam_theta_mom_xB_${xBbin}_Q2_${Q2bin}", h_theta_mom).fill(gam.p(), Math.toDegrees(gam.theta()))
          hists.computeIfAbsent("/dvcs/corr/gam_theta_phi_xB_${xBbin}_Q2_${Q2bin}", h_theta_phi).fill(gam_phi_convention, Math.toDegrees(gam.theta()))

          hists.computeIfAbsent("/dvcs/corr/prot_theta_t_xB_${xBbin}_Q2_${Q2bin}", h_theta_t).fill(t, Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/dvcs/corr/prot_phi_t_xB_${xBbin}_Q2_${Q2bin}", h_phi_t).fill(t, pro_phi_convention)
          hists.computeIfAbsent("/dvcs/corr/prot_theta_trento_xB_${xBbin}_Q2_${Q2bin}", h_theta_trento).fill(TrentoAng, Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/dvcs/corr/prot_phi_trento_xB_${xBbin}_Q2_${Q2bin}", h_phi_trento).fill(TrentoAng, pro_phi_convention)
          hists.computeIfAbsent("/dvcs/corr/gam_theta_t_xB_${xBbin}_Q2_${Q2bin}", h_theta_t).fill(t, Math.toDegrees(gam.theta()))
          hists.computeIfAbsent("/dvcs/corr/gam_phi_t_xB_${xBbin}_Q2_${Q2bin}", h_phi_t).fill(t, gam_phi_convention)
          hists.computeIfAbsent("/dvcs/corr/gam_theta_trento_xB_${xBbin}_Q2_${Q2bin}", h_theta_trento).fill(TrentoAng, Math.toDegrees(gam.theta()))
          hists.computeIfAbsent("/dvcs/corr/gam_phi_trento_xB_${xBbin}_Q2_${Q2bin}", h_phi_trento).fill(TrentoAng, gam_phi_convention)


          hists.computeIfAbsent("/dvcs/elec_polar_sec"+ele_sec, h_polar_rate).fill(Math.toDegrees(ele.theta()))
          hists.computeIfAbsent("/dvcs/prot_polar", h_polar_rate).fill(Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/dvcs/gam_polar", h_polar_rate).fill(Math.toDegrees(gam.theta()))


          hists.computeIfAbsent("/dvcs/tdep/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
          hists.computeIfAbsent("/events/events", h_events).fill(3.5)  
          if (event.status[dsets.pindex[1]]<4000 && event.status[dsets.pindex[1]]>=2000){
            hists.computeIfAbsent("/dvcs/tdep/gam_fd/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
          hists.computeIfAbsent("/events/events", h_events).fill(4.5)  
          }
          else if (event.status[dsets.pindex[1]]>=1000){
            hists.computeIfAbsent("/dvcs/tdep/gam_ft/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
            hists.computeIfAbsent("/events/events", h_events).fill(5.5)  
          }

          if (event.status[dsets.pindex[1]]>=4000){
            hists.computeIfAbsent("/dvcs/tdep/pro_cd/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
            hists.computeIfAbsent("/events/events", h_events).fill(6.5)  
            if (event.status[dsets.pindex[1]]<4000 && event.status[dsets.pindex[1]]>=2000){
              hists.computeIfAbsent("/dvcs/tdep/pro_cd/gam_fd/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
              hists.computeIfAbsent("/events/events", h_events).fill(7.5)  
            }
            else if (event.status[dsets.pindex[1]]>=1000){
              hists.computeIfAbsent("/dvcs/tdep/pro_cd/gam_ft/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
              hists.computeIfAbsent("/events/events", h_events).fill(8.5)  
            }
          }
          else if (event.status[dsets.pindex[1]]>=2000){
            hists.computeIfAbsent("/dvcs/tdep/pro_fd/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
            hists.computeIfAbsent("/events/events", h_events).fill(9.5)  
            if (event.status[dsets.pindex[1]]<4000 && event.status[dsets.pindex[1]]>=2000){
              hists.computeIfAbsent("/dvcs/tdep/pro_fd/gam_fd/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
              hists.computeIfAbsent("/events/events", h_events).fill(10.5)  
            }
            else if (event.status[dsets.pindex[1]]>=1000){
              hists.computeIfAbsent("/dvcs/tdep/pro_fd/gam_ft/h_xB_${xBbin}_Q2_${Q2bin}", h_t).fill(t)
              hists.computeIfAbsent("/events/events", h_events).fill(11.5)  
            }
          }

          // fill phi dependence on 3 fold binning (xB, Q2, t)
          def bin_number = binnumber(xB, ele.theta(), t)
          hists.computeIfAbsent("/dvcs/h_phi_bin_$bin_number", h_cross_section).fill(TrentoAng)
          hists.computeIfAbsent("/dvcs/h_Q2_xB_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
          
          if (event.status[dsets.pindex[1]]>=4000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_pro_CD_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_pro_CD_bin_$bin_number", h_cross_section).fill(TrentoAng)
          }

          if (event.status[dsets.pindex[2]]<2000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_gam_FT_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_gam_FT_bin_$bin_number", h_cross_section).fill(TrentoAng)
          }

          if (event.status[dsets.pindex[1]]>=4000 && event.status[dsets.pindex[2]]<2000){
            hists.computeIfAbsent("/dvcs/h_Q2_xB_pro_CD_gam_FT_bin_$bin_number", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/h_phi_pro_CD_gam_FT_bin_$bin_number", h_cross_section).fill(TrentoAng)
            hists.computeIfAbsent("/dvcs/hmm2_epg", hmm2_2).fill(VMISS.mass2())
          }            

          if (event.status[dsets.pindex[1]]>=4000){
            hists.computeIfAbsent("/dvcs/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/dvcs/prot_azimuth_CD", h_azimuth_rate).fill(pro_phi)
          }
          else if (event.status[dsets.pindex[1]]<4000){
            hists.computeIfAbsent("/dvcs/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/dvcs/prot_azimuth_FD", h_azimuth_rate).fill(pro_phi)
          }
        } // exclusivity cuts ended
        //add here for analysis
      }// event with e, p, g
    }// event with particles
  }
}