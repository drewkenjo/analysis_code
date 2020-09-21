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

class dvcs_EB{

  //defining histograms
  def hists = new ConcurrentHashMap()

  // missing mass
  def h_mm2 = {new H1F("$it", "$it", 1500, -2, 4)}
  def h_mm2_2 = {new H1F("$it", "$it", 100, -0.2, 0.2)}
  def h_mm2_me = {new H2F("$it", "$it", 1000, -1, 3, 100, -0.1, 0.1)}

  // invaraiant mass
  def h_inv_mass_gg = {new H1F("$it", "$it", 1000, 0, 0.2)}

  // angle between planes
  def h_angle = {new H1F("$it", "$it", 1900, -5 ,185)}

  // kinematic variables (correlation)
  def h_theta_mom = {new H2F("$it", "$it", 120, 0, 12, 100, 0, 100)}
  def h_phi_mom = {new H2F("$it", "$it", 120, 0, 12, 360, -180, 180)}
  def h_theta_phi = {new H2F("$it", "$it", 360, -180, 180, 100, 0, 100)}
  def h_theta_t = {new H2F("$it", "$it", 100, 0, 4, 100, 0, 100)}
  def h_phi_t = {new H2F("$it", "$it", 100, 0, 4, 360, -180, 180)}
  def h_theta_trento = {new H2F("$it", "$it", 360, 0, 360, 100, 0, 100)}
  def h_phi_trento = {new H2F("$it", "$it", 360, 0, 360, 360, -180, 180)}

  def h_t_trento =  {new H2F("$it", "$it", 360, 0, 360, 100, 0 , 4)}
  def h_t_trento_logarithmic =  {new H2F("$it", "$it", 360, 0, 360, 30, -4 , 0.6)}

  def h_t_t = {new H2F("$it", "$it", 100,0,4, 100,0,4)}

  def h_Q2_xB_logarithmic = {new H2F("$it", "$it", 30, -1.31, -0.096, 30, 0, 1)}
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
  def h_Q2_W = {new H2F("$it","$it",100, 0, 10, 100, 0, 12)}

  //binning
  def xB_array = [0, 0.1, 0.15, 0.2, 0.3, 0.4, 1]
  def t_array = [0, 0.09, 0.13, 0.18, 0.23, 0.3, 0.39, 0.52, 0.72, 1.1, 2]
  def Q2_array = [1, 1.7, 2.22, 2.7, 3.5, 11]

  def xB_bin = {xB ->
    int xBbin = xB_array.findIndexOf{ xB < it} -1
    if (xBbin == -2) xBbin = xBbin + xB_array.size() +1 //length
    return xBbin
  }

  def Q2_bin = {Q2 ->
    int Q2bin = Q2_array.findIndexOf{ Q2 < it} -1
    if (Q2bin == -2) Q2bin = Q2bin + Q2_array.size() +1 //length
    return Q2bin
  }

  def t_bin = {t ->
    int tbin = t_array.findIndexOf{ t < it} -1
    if (tbin == -2) tbin = tbin + t_array.size() +1 //length
    return tbin
  }

  // pid histograms
  def h_vz_mom = {new H2F("$it", "$it",110,0,11,800,-40,40)}     //vz vs p
  def h_vz_theta = {new H2F("$it", "$it",100,0,100,800,-40,40)} //vz vs p
  def h_pcalecal = {new H2F("$it","$it",100,0,0.5,100,0,0.5)}   //ecal/p vs pcal/p
  def h_Sampl_mom =  {new H2F("$it", "$it",110,0,11,100,0,1)}   //Sampling fraction vs p

  def h_theta_phi_ft = {new H2F("$it", "$it", 360, -180, 180, 100, 1, 6)}
  def h_theta_mom_ft = {new H2F("$it", "$it", 120, 0, 12, 100, 1, 6)}

  def phi_convention = {phi ->
    phi=phi + 180
    if (phi>180) phi = phi - 360
    if (phi<-180) phi = phi + 360
    return phi
  }

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
      def dsets = DVCS.getEPG_EB(event)
      def (ele, pro, gam) = dsets*.particle.collect{it ? it.vector() : null} 
      // process only if there's a epg set in coincidence
      if(ele!=null) {
        // event number histograms
        hists.computeIfAbsent("/events/events", h_events).fill(2.5)  

        // get pindex, sector, status
        def (ele_ind, pro_ind, gam_ind) = dsets*.pindex
        def (ele_sec, pro_sec, gam_sec) = dsets*.sector
        def (ele_status, pro_status, gam_status) = dsets*.status
        
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
        def nu = KinTool.calcNu(beam, ele)
        def M = PDGDatabase.getParticleMass(2212)
        def costheta = VGS.vect().dot(gam.vect())/VGS.vect().mag()/gam.vect().mag()
        def t2 = (M*Q2+2*M*nu*(nu-Math.sqrt(nu*nu+Q2)*costheta))/(M+nu-Math.sqrt(nu*nu+Q2)*costheta)

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
        hists.computeIfAbsent("/epg/h_mm2_ep", h_mm2).fill(VmissG.mass2())
        hists.computeIfAbsent("/epg/h_mm2_eg", h_mm2).fill(VmissP.mass2())
        hists.computeIfAbsent("/epg/h_mm2_epg", h_mm2_2).fill(VMISS.mass2())

        // kinematic range
        hists.computeIfAbsent("/epg/h_Q2_xB", h_Q2_xB).fill(xB,Q2)
        hists.computeIfAbsent("/epg/h_Q2_t", h_Q2_t).fill(t,Q2)
        hists.computeIfAbsent("/epg/h_t_xB", h_t_xB).fill(xB,t)
        hists.computeIfAbsent("/epg/h_Q2_theta", h_Q2_theta).fill(Math.toDegrees(ele.theta()),Q2);
        hists.computeIfAbsent("/epg/h_Q2_xB_sec"+ele_sec, h_Q2_xB).fill(xB,Q2)
        hists.computeIfAbsent("/epg/h_W_sec"+ele_sec, h_W).fill(W)
        hists.computeIfAbsent("/epg/h_t_sec"+ele_sec, h_t).fill(t)
        hists.computeIfAbsent("/epg/h_phi_sec"+ele_sec, h_cross_section).fill(TrentoAng) 
        hists.computeIfAbsent("/epg/h_y_sec"+ele_sec, h_y).fill(KinTool.calcY(beam, ele))
        hists.computeIfAbsent("/epg/h_Q2_W_sec"+ele_sec, h_Q2_W).fill(W, Q2)

        // working on binning // don't have to fill this because binning is over.
        // double theta_e = Math.toDegrees(ele.theta()) /5
        // String theta_label = String.format("%.0f_theta_%.0f",5*theta_e.trunc(),5*theta_e.trunc()+5)
        // h_Q2_xB_cond[theta_label].fill(xB,Q2)
        if (W>2)  hists.computeIfAbsent("/epg/h_Q2_xB_W>2", h_Q2_xB).fill(xB,Q2)
        else hists.computeIfAbsent("/epg/h_Q2_xB_W<2", h_Q2_xB).fill(xB,Q2)

        // check CD alignment
        if (pro_status>=4000){
          hists.computeIfAbsent("/epg/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/epg/prot_azimuth_CD", h_azimuth_rate).fill(pro_phi)
          if (W>2){
            hists.computeIfAbsent("/epg/prot_polar_CD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/epg/prot_azimuth_CD_W>2", h_azimuth_rate).fill(pro_phi)
          }
        }
        else if (pro_status<4000){
          hists.computeIfAbsent("/epg/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
          hists.computeIfAbsent("/epg/prot_azimuth_FD", h_azimuth_rate).fill(pro_phi)
          if (W>2){
            hists.computeIfAbsent("/epg/prot_polar_FD_W>2", h_polar_rate).fill(Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/epg/prot_azimuth_FD_W>2", h_azimuth_rate).fill(pro_phi)
          }
        }

        def eidep=0
        def eodep = 0
        def pcaldep = 0
        if( event.ecal_inner_status.contains(ele_ind) )  eidep = event.ecal_inner_energy[ele_ind]
        if( event.ecal_outer_status.contains(ele_ind) )  eodep = event.ecal_outer_energy[ele_ind]
        if( event.pcal_status.contains(ele_ind) )        pcaldep = event.pcal_energy[ele_ind]
        def edep = eidep + eodep + pcaldep

        //electron pid
        hists.computeIfAbsent("/epg/pid/ele/vz_mom_S"+ele_sec, h_vz_mom).fill(ele.p(), event.vz[ele_ind])
        hists.computeIfAbsent("/epg/pid/ele/vz_theta_S"+ele_sec, h_vz_theta).fill(Math.toDegrees(ele.theta()), event.vz[ele_ind])
        hists.computeIfAbsent("/epg/pid/ele/pcalecal_S"+ele_sec, h_pcalecal).fill(eidep/ele.p(), pcaldep/ele.p())
        hists.computeIfAbsent("/epg/pid/ele/Sampl_mom_S"+ele_sec, h_Sampl_mom).fill(ele.p(), edep/ele.p())
        hists.computeIfAbsent("/epg/pid/ele/theta_phi_S"+ele_sec, h_theta_phi).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))
        hists.computeIfAbsent("/epg/pid/ele/theta_mom_S"+ele_sec, h_theta_mom).fill(ele.p(), Math.toDegrees(ele.theta()))

        //proton pid
        def pro_string = ''
        if (pro_status>=4000) pro_string = 'cd'
        else if (pro_status>=2000) pro_string = 'fd_S'+pro_sec 

        hists.computeIfAbsent("/epg/pid/pro/vz_mom_"+pro_string, h_vz_mom).fill(pro.p(), event.vz[pro_ind])
        hists.computeIfAbsent("/epg/pid/pro/vzdiff_mom_"+pro_string, h_vz_mom).fill(pro.p(), event.vz[pro_ind]-event.vz[ele_ind])
        hists.computeIfAbsent("/epg/pid/pro/vz_theta_"+pro_string, h_vz_theta).fill(Math.toDegrees(pro.theta()), event.vz[pro_ind])
        hists.computeIfAbsent("/epg/pid/pro/theta_phi_"+pro_string, h_theta_phi).fill(Math.toDegrees(pro.phi()), Math.toDegrees(pro.theta()))
        hists.computeIfAbsent("/epg/pid/pro/theta_mom_"+pro_string, h_theta_mom).fill(pro.p(), Math.toDegrees(pro.theta()))

        //gamma pid
        def gam_string = ''
        if (gam_status<4000 && gam_status>=2000) gam_string = 'fd_S'+gam_sec
        else if (gam_status>=1000 && gam_status<2000) gam_string = 'ft' 

        if (gam_string == 'ft'){
          hists.computeIfAbsent("/epg/pid/gam/theta_phi_"+gam_string, h_theta_phi_ft).fill(Math.toDegrees(gam.phi()), Math.toDegrees(gam.theta()))
          hists.computeIfAbsent("/epg/pid/gam/theta_mom_"+gam_string, h_theta_mom_ft).fill(gam.p(), Math.toDegrees(gam.theta()))
        }
        else{
          hists.computeIfAbsent("/epg/pid/gam/theta_phi_"+gam_string, h_theta_phi).fill(Math.toDegrees(gam.phi()), Math.toDegrees(gam.theta()))
          hists.computeIfAbsent("/epg/pid/gam/theta_mom_"+gam_string, h_theta_mom).fill(gam.p(), Math.toDegrees(gam.theta()))
        }

        if (DVCS.KineCuts(xB, Q2, W, ele, gam)){

          hists.computeIfAbsent("/events/events", h_events).fill(3.5)  

          //excl cuts
          hists.computeIfAbsent("/excl/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
          hists.computeIfAbsent("/excl/missing_energy", h_mm2).fill(VMISS.e())
          hists.computeIfAbsent("/excl/missing_mass_epg", h_mm2).fill(VMISS.mass2())
          hists.computeIfAbsent("/excl/missing_mass_eg", h_mm2).fill(VmissP.mass2())
          hists.computeIfAbsent("/excl/missing_mass_ep", h_mm2).fill(VmissG.mass2())
          hists.computeIfAbsent("/excl/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
          hists.computeIfAbsent("/excl/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
          hists.computeIfAbsent("/excl/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
          hists.computeIfAbsent("/excl/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())

          if (KinTool.Vangle(gam.vect(),ele.vect())>4){
            hists.computeIfAbsent("/excl/cuts/cone_angle/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/cone_angle/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/cone_angle/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/cone_angle/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/cone_angle/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/cone_angle/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/cone_angle/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/cone_angle/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/cone_angle/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (VMISS.e()<1.5){
            hists.computeIfAbsent("/excl/cuts/missing_energy/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_energy/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/missing_energy/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_energy/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_energy/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_energy/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/missing_energy/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_energy/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/missing_energy/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (VMISS.mass2() <0.2 && VMISS.mass2() >-0.2 ){
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/missing_mass_epg/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (VmissP.mass2() < 3 && VmissP.mass2() > -0.25){
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/missing_mass_eg/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (VmissG.mass2() < 1 && VmissG.mass2() > -1){
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/missing_mass_ep/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 0.3){
            hists.computeIfAbsent("/excl/cuts/missing_pt/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_pt/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/missing_pt/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_pt/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_pt/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/missing_pt/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/missing_pt/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/missing_pt/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/missing_pt/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (KinTool.Vangle(gam.vect(),VmissG.vect()) < 3){
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/recon_gam_cone_angle/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
          }

          if (KinTool.Vangle(Vhad2,Vhadr) < 4){
            hists.computeIfAbsent("/excl/cuts/coplanarity/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/coplanarity/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/coplanarity/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/coplanarity/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/coplanarity/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/coplanarity/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/coplanarity/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/coplanarity/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/coplanarity/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())
         }

          if (DVCS.ExclCuts(gam, ele, VMISS, VmissP, VmissG, Vhadr, Vhad2)){

            hists.computeIfAbsent("/events/events", h_events).fill(4.5)  
            
            //calc tcol tmin
            def E = 10.6
            def tmin = M*M*xB*xB/(1-xB+xB*M*M/Q2)
            def tcol = Q2*(Q2-2*xB*M*E)/xB/(Q2-2*M*E)
            // fill t dependence on 2 fold binning (xB, Q2)
            int xBbin = 1 + 2 * Math.floor(xB/0.2)
            int Q2bin = 1 + 2 * Math.floor(Q2/2)

            //electron pid
            hists.computeIfAbsent("/excl/pid/ele/vz_mom_S"+ele_sec, h_vz_mom).fill(ele.p(), event.vz[ele_ind])
            hists.computeIfAbsent("/excl/pid/ele/vz_theta_S"+ele_sec, h_vz_theta).fill(Math.toDegrees(ele.theta()), event.vz[ele_ind])
            hists.computeIfAbsent("/excl/pid/ele/pcalecal_S"+ele_sec, h_pcalecal).fill(eidep/ele.p(), pcaldep/ele.p())
            hists.computeIfAbsent("/excl/pid/ele/Sampl_mom_S"+ele_sec, h_Sampl_mom).fill(ele.p(), edep/ele.p())
            hists.computeIfAbsent("/excl/pid/ele/theta_phi_S"+ele_sec, h_theta_phi).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))
            hists.computeIfAbsent("/excl/pid/ele/theta_mom_S"+ele_sec, h_theta_mom).fill(ele.p(), Math.toDegrees(ele.theta()))

            //proton pid
            hists.computeIfAbsent("/excl/pid/pro/vz_mom_"+pro_string, h_vz_mom).fill(pro.p(), event.vz[pro_ind])
            hists.computeIfAbsent("/excl/pid/pro/vzdiff_mom_"+pro_string, h_vz_mom).fill(pro.p(), event.vz[pro_ind]-event.vz[ele_ind])
            hists.computeIfAbsent("/excl/pid/pro/vz_theta_"+pro_string, h_vz_theta).fill(Math.toDegrees(pro.theta()), event.vz[pro_ind])
            hists.computeIfAbsent("/excl/pid/pro/theta_phi_"+pro_string, h_theta_phi).fill(Math.toDegrees(pro.phi()), Math.toDegrees(pro.theta()))
            hists.computeIfAbsent("/excl/pid/pro/theta_mom_"+pro_string, h_theta_mom).fill(pro.p(), Math.toDegrees(pro.theta()))

            //gamma pid
            if (gam_string == 'ft'){
              hists.computeIfAbsent("/excl/pid/gam/theta_phi_"+gam_string, h_theta_phi_ft).fill(Math.toDegrees(gam.phi()), Math.toDegrees(gam.theta()))
              hists.computeIfAbsent("/excl/pid/gam/theta_mom_"+gam_string, h_theta_mom_ft).fill(gam.p(), Math.toDegrees(gam.theta()))
            }
            else{
              hists.computeIfAbsent("/excl/pid/gam/theta_phi_"+gam_string, h_theta_phi).fill(Math.toDegrees(gam.phi()), Math.toDegrees(gam.theta()))
              hists.computeIfAbsent("/excl/pid/gam/theta_mom_"+gam_string, h_theta_mom).fill(gam.p(), Math.toDegrees(gam.theta()))
            }

            //excl cuts
            hists.computeIfAbsent("/excl/cuts/all/cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),ele.vect()))
            hists.computeIfAbsent("/excl/cuts/all/missing_energy", h_mm2).fill(VMISS.e())
            hists.computeIfAbsent("/excl/cuts/all/missing_mass_epg", h_mm2).fill(VMISS.mass2())
            hists.computeIfAbsent("/excl/cuts/all/missing_mass_eg", h_mm2).fill(VmissP.mass2())
            hists.computeIfAbsent("/excl/cuts/all/missing_mass_ep", h_mm2).fill(VmissG.mass2())
            hists.computeIfAbsent("/excl/cuts/all/missing_pt", h_mm2).fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()))
            hists.computeIfAbsent("/excl/cuts/all/recon_gam_cone_angle", h_angle).fill(KinTool.Vangle(gam.vect(),VmissG.vect()))
            hists.computeIfAbsent("/excl/cuts/all/coplanarity", h_angle).fill(KinTool.Vangle(Vhad2,Vhadr))
            hists.computeIfAbsent("/excl/cuts/all/mm2epg_me", h_mm2_me).fill(VMISS.e(), VMISS.mass2())

            hists.computeIfAbsent("/dvcs/corr/tmin", h_Q2_xB).fill(xB,Q2,tmin)
            hists.computeIfAbsent("/dvcs/corr/tcol", h_Q2_xB).fill(xB,Q2,tcol)
            hists.computeIfAbsent("/dvcs/binning/h_Q2_xB", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/binning/h_t_trento", h_t_trento).fill(TrentoAng, t)
            hists.computeIfAbsent("/dvcs/binning/h_Q2_t", h_Q2_t).fill(t,Q2)
            hists.computeIfAbsent("/dvcs/binning/h_Q2_xB_logarithmic", h_Q2_xB_logarithmic).fill(Math.log10(xB), Math.log10(Q2), 1/Q2/xB)
            hists.computeIfAbsent("/dvcs/binning/h_t_trento_logarithmic", h_t_trento_logarithmic).fill(TrentoAng, Math.log10(t2), 1/t2)

            hists.computeIfAbsent("/dvcs/binning/h_t_xB", h_t_xB).fill(xB,t)
            hists.computeIfAbsent("/dvcs/binning/h_Q2_theta", h_Q2_theta).fill(Math.toDegrees(ele.theta()),Q2);
            hists.computeIfAbsent("/dvcs/binning/h_Q2_xB_sec"+ele_sec, h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/binning/h_W_sec"+ele_sec, h_W).fill(W)
            hists.computeIfAbsent("/dvcs/binning/h_t_sec"+ele_sec, h_t).fill(t)
            hists.computeIfAbsent("/dvcs/binning/h_phi_sec"+ele_sec, h_cross_section).fill(TrentoAng) 
            hists.computeIfAbsent("/dvcs/binning/h_y_sec"+ele_sec, h_y).fill(KinTool.calcY(beam, ele))
            hists.computeIfAbsent("/dvcs/binning/h_Q2_W_sec"+ele_sec, h_Q2_W).fill(W, Q2)

            def pro_phi_convention = phi_convention(Math.toDegrees(pro.phi()-ele.phi()))
            def gam_phi_convention = phi_convention(Math.toDegrees(gam.phi()-ele.phi()))

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

            hists.computeIfAbsent("/dvcs/corr/h_t_t", h_t_t).fill(t, t2)
            hists.computeIfAbsent("/dvcs/corr/h_Q2_xB_t1", h_Q2_xB).fill(xB,Q2,t)
            hists.computeIfAbsent("/dvcs/corr/h_Q2_xB_t2", h_Q2_xB).fill(xB,Q2,t2)

            // check CD alignment
            if (pro_status>=4000){
              hists.computeIfAbsent("/events/events", h_events).fill(5.5)  

              hists.computeIfAbsent("/dvcs/prot_polar_CD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
              hists.computeIfAbsent("/dvcs/prot_azimuth_CD", h_azimuth_rate).fill(pro_phi)
            }
            else if (pro_status<4000 && pro_status>2000){
              hists.computeIfAbsent("/events/events", h_events).fill(6.5)
                
              hists.computeIfAbsent("/dvcs/prot_polar_FD", h_polar_rate).fill(Math.toDegrees(pro.theta()))
              hists.computeIfAbsent("/dvcs/prot_azimuth_FD", h_azimuth_rate).fill(pro_phi)
            }

            def number_of_photons = gamma_selector.applyCuts_Stefan(event).size()
            hists.computeIfAbsent("/dvcs/number_of_photons", h_events).fill(number_of_photons)
            if (number_of_photons>1){
              def ind_gam2 = gamma_selector.applyCuts_Stefan(event).max{ind->
                if (ind!=dsets.pindex[2]) new Vector3(*[event.px, event.py, event.pz].collect{it[ind]}).mag2()}
              def gam2 = LorentzVector.withPID(22, *[event.px, event.py, event.pz].collect{it[ind_gam2]})
              hists.computeIfAbsent("/dvcs/pi0/h_inv_mass_gg", h_inv_mass_gg).fill((gam + gam2).mass())
            }

            def xBbin2 = xB_bin(xB)
            def Q2bin2 = Q2_bin(Q2)
            def tbin = t_bin(t2)
            def helicity = event.helicity

            hists.computeIfAbsent("/dvcs/heli_$helicity/h_Q2_xB_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_Q2_xB).fill(xB,Q2)
            hists.computeIfAbsent("/dvcs/heli_$helicity/h_trento_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_cross_section).fill(TrentoAng)
            
            if (pro_status>=4000){
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_Q2_xB_pro_CD_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_Q2_xB).fill(xB,Q2)
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_trento_pro_CD_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_cross_section).fill(TrentoAng)
            }

            if (gam_status<2000){
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_Q2_xB_gam_FT_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_Q2_xB).fill(xB,Q2)
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_trento_gam_FT_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_cross_section).fill(TrentoAng)
            }

            if (pro_status>=4000 && gam_status<2000){
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_Q2_xB_pro_CD_gam_FT_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_Q2_xB).fill(xB,Q2)
              hists.computeIfAbsent("/dvcs/heli_$helicity/h_trento_pro_CD_gam_FT_xB_${xBbin2}_Q2_${Q2bin2}_t_${tbin}", h_cross_section).fill(TrentoAng)
            }            
          } // exclusivity cuts ended
        }//kine cuts ended
      }// event with e, p, g
    }// event with particles
  }
}