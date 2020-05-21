package pid.gamma

import org.jlab.detector.base.DetectorType
import event.Event

class GammaFromEvent {

  def ebeam = 10.6
  def ebPID = 22

  // Stefan's original comments below
  /// 1.2 PCAL fiducial cut based on fitted sampling fraction (I woudl recommend this cut for photons)
  ///     For this cut, the cut criterium is teh drop of the sampling fraction 
  ///     loose: The mean value of E/p is allowed to drop to 2 RMS of the distribution 
  ///     medium: The mean value of E/p is allowed to drop to 1.5 RMS of the distribution 
  ///     tight: The mean value of E/p is allowed to drop to 1 RMS of the distribution 

  // 2 sigma (inb adjusted):

  def par_0_min, par_0_min_inb, par_0_min_outb 
  def par_1_min, par_1_min_inb, par_1_min_outb 
  def par_2_min, par_2_min_inb, par_2_min_outb
  def par_3_min, par_3_min_inb, par_3_min_outb 
  def par_4_min, par_4_min_inb, par_4_min_outb 
  def par_0_max, par_0_max_inb, par_0_max_outb 
  def par_1_max, par_1_max_inb, par_1_max_outb 
  def par_2_max, par_2_max_inb, par_2_max_outb 
  def par_3_max, par_3_max_inb, par_3_max_outb 
  def par_4_max, par_4_max_inb, par_4_max_outb 

  def par_0_min_2sigma_inb = [13.771, 25.639, 28.4616, 34.2333, 41.777, 18.041]
  def par_1_min_2sigma_inb = [-14.3952, -27.4437, -29.5074, -40.9268, -33.056, -19.1539]
  def par_2_min_2sigma_inb = [-0.0579567, 2.12184, 2.7033, 4.72287, 1.80765, 0.648548]
  def par_3_min_2sigma_inb = [0.0190146, -0.0315046, -0.0523381, -0.0919215, -0.00542624, 0.00836905]
  def par_4_min_2sigma_inb = [-0.000222315, 0.000243574, 0.000453275, 0.000771317, -0.000150837, -0.000208545]
  def par_0_max_2sigma_inb = [-19.2009, -16.4848, -47.8295, -24.0029, -25.096, -19.2967]
  def par_1_max_2sigma_inb = [19.3148, 11.0556, 47.4188, 23.6525, 19.3032, 19.9627]
  def par_2_max_2sigma_inb = [-0.66582, 1.35067, -5.54184, -1.20742, 0.0744728, -0.714339]
  def par_3_max_2sigma_inb = [-0.00592537, -0.0623736, 0.123022, 0.00276304, -0.0334298, -0.0077081]
  def par_4_max_2sigma_inb = [0.000187643, 0.000759023, -0.00120291, 0.000128345, 0.000486216, 0.000224336]

  def par_0_min_2sigma_outb = [64.895, 69.4726, 37.4087, 57.2897, 36.5138, 66.2482]
  def par_1_min_2sigma_outb = [-49.7813, -52.7634, -28.7868, -44.398, -27.2841, -50.6646]
  def par_2_min_2sigma_outb = [3.79889, 4.06934, 1.52759, 3.26992, 1.25399, 3.8817]
  def par_3_min_2sigma_outb = [-0.0389919, -0.0418169, -0.0120817, -0.0329207, -0.00772183, -0.0398076]
  def par_4_min_2sigma_outb = [0,0,0,0,0,0]
  def par_0_max_2sigma_outb = [-58.8252, -42.022, -35.6843, -62.0889, -25.7336, -62.4078]
  def par_1_max_2sigma_outb = [46.3788, 32.7105, 29.2649, 48.8274, 21.4091, 49.5489]
  def par_2_max_2sigma_outb = [-3.47241, -2.00905, -1.83462, -3.7453, -0.830711, -3.86809]
  def par_3_max_2sigma_outb = [0.0350037, 0.018892, 0.0188074, 0.0381712, 0.00422423, 0.040186]
  def par_4_max_2sigma_outb = [0,0,0,0,0,0]

  // 1.5 sigma (inb adjusted, 4 min replaced by 2 min!):

  def par_0_min_15sigma_inb = [25.2996, 19.3705, 59.5003, 19.3705, 26.9823, 21.8217]
  def par_1_min_15sigma_inb = [-26.1158, -19.5271, -55.9639, -19.5271, -22.4489, -23.0262]
  def par_2_min_15sigma_inb = [2.09145, 0.72118, 6.4372, 0.72118, 0.890449, 1.48469]
  def par_3_min_15sigma_inb = [-0.041483, 0.00293465, -0.13059, 0.00293465, -0.000774943, -0.0190315]
  def par_4_min_15sigma_inb = [0.000456261, -0.000109323, 0.00115246, -0.000109323, -7.12074e-05, 0.000125463]
  def par_0_max_15sigma_inb = [-26.7425, -15.9009, -47.556, -21.5038, -33.9197, -24.0325]
  def par_1_max_15sigma_inb = [26.004, 10.5989, 41.6295, 21.1734, 31.1811, 23.3122]
  def par_2_max_15sigma_inb = [-1.76638, 1.10168, -3.66934, -0.969572, -2.51229, -1.23308]
  def par_3_max_15sigma_inb = [0.0227414, -0.0455311, 0.0493171, 0.00373945, 0.0443308, 0.00762314]
  def par_4_max_15sigma_inb = [-0.000111102, 0.000503536, -0.000151053, 5.22425e-05, -0.000403627, 4.24553e-05]

  def par_0_min_15sigma_outb = [57.0314, 70.411, 74.9683, 53.9146, 44.7614, 64.012]
  def par_1_min_15sigma_outb = [-43.0803, -53.3163, -59.0842, -41.6436, -35.2193, -48.8726]
  def par_2_min_15sigma_outb = [2.99452, 4.11397, 5.27234, 2.96164, 2.44681, 3.67536]
  def par_3_min_15sigma_outb = [-0.0287862, -0.042257, -0.0638554, -0.0293148, -0.0264986, -0.0372214]
  def par_4_min_15sigma_outb = [0,0,0,0,0,0]
  def par_0_max_15sigma_outb = [-52.0537, -48.3703, -94.6197, -54.6123, -79.9164, -55.3222]
  def par_1_max_15sigma_outb = [40.7573, 38.333, 73.5425, 42.9251, 64.9277, 42.8186]
  def par_2_max_15sigma_outb = [-2.8105, -2.83403, -6.88649, -3.08431, -6.18973, -2.99337]
  def par_3_max_15sigma_outb = [0.0266371, 0.0318955, 0.0851474, 0.0301575, 0.0787016, 0.0286132]
  def par_4_max_15sigma_outb = [0,0,0,0,0,0]

  // 1 sigma (inb adjusted):

  def par_0_min_1sigma_inb = [34.1128, 26.6064, 65.3241, 44.0344, 54.5539, 25.7356]
  def par_1_min_1sigma_inb = [-30.7179, -26.3373, -58.7761, -35.918, -47.3194, -25.3968]
  def par_2_min_1sigma_inb = [2.31272, 1.85141, 6.48495, 2.34733, 4.58872, 1.76128]
  def par_3_min_1sigma_inb = [-0.0351384, -0.023687, -0.121042, -0.0170119, -0.0778135, -0.0243688]
  def par_4_min_1sigma_inb = [0.000262438, 0.000120765, 0.000936822, -7.66933e-05, 0.000539922, 0.000156061]
  def par_0_max_1sigma_inb = [-31.6364, -28.7094, -35.6017, -30.1334, -61.5491, -30.9496]
  def par_1_max_1sigma_inb = [27.253, 20.8471, 26.4139, 27.7852, 55.5266, 28.8408]
  def par_2_max_1sigma_inb = [-1.53381, -0.254236, -0.77312, -1.85849, -6.03473, -2.07467]
  def par_3_max_1sigma_inb = [0.00817262, -0.0259995, -0.0300504, 0.0211806, 0.113112, 0.0285026]
  def par_4_max_1sigma_inb = [0.000129402, 0.000495665, 0.000740296, -6.30543e-05, -0.000871381, -0.000162669]

  def par_0_min_1sigma_outb = [58.8694, 75.494, 119.951, 52.9731, 111.593, 53.6272]
  def par_1_min_1sigma_outb = [-43.6605, -55.9065, -89.9595, -41.2703, -85.093, -40.5409]
  def par_2_min_1sigma_outb = [3.03906, 4.28761, 8.38977, 3.00994, 8.06472, 2.75155]
  def par_3_min_1sigma_outb = [-0.0301444, -0.0440019, -0.0999269, -0.0313006, -0.0984858, -0.0264477]
  def par_4_min_1sigma_outb = [0,0,0,0,0,0]
  def par_0_max_1sigma_outb = [-40.256, -58.3938, -60.3614, -57.7244, -102.98, -51.0424]
  def par_1_max_1sigma_outb = [31.4367, 43.3923, 45.8203, 42.9619, 79.613, 38.5613]
  def par_2_max_1sigma_outb = [-1.78797, -3.14225, -3.49161, -2.9105, -7.51346, -2.46405]
  def par_3_max_1sigma_outb = [0.0146775, 0.0334355, 0.0366689, 0.0277125, 0.0916122, 0.0223185]
  def par_4_max_1sigma_outb = [0,0,0,0,0,0]


  def gam_cut_strictness_lvl = null
  def pcal_fiducial_cut_lvl = null

  void setGammaCutStrictness(gam_cut_strictness){
    gam_cut_strictness_lvl=gam_cut_strictness
    pcal_fiducial_cut_lvl=gam_cut_strictness["pcal_fiducial_cut_lvl"] 

    println("[GammaFromEvent::setGammaCutStrictness] -> pcal_fiducial_cut_lvl " + pcal_fiducial_cut_lvl)

    //loose
    if (pcal_fiducial_cut_lvl==0){
      par_0_min_inb = par_0_min_1sigma_inb 
      par_1_min_inb = par_1_min_1sigma_inb 
      par_2_min_inb = par_2_min_1sigma_inb
      par_3_min_inb = par_3_min_1sigma_inb 
      par_4_min_inb = par_4_min_1sigma_inb 
      par_0_max_inb = par_0_max_1sigma_inb 
      par_1_max_inb = par_1_max_1sigma_inb 
      par_2_max_inb = par_2_max_1sigma_inb 
      par_3_max_inb = par_3_max_1sigma_inb 
      par_4_max_inb = par_4_max_1sigma_inb 

      par_0_min_outb = par_0_min_1sigma_outb 
      par_1_min_outb = par_1_min_1sigma_outb 
      par_2_min_outb = par_2_min_1sigma_outb
      par_3_min_outb = par_3_min_1sigma_outb 
      par_4_min_outb = par_4_min_1sigma_outb 
      par_0_max_outb = par_0_max_1sigma_outb 
      par_1_max_outb = par_1_max_1sigma_outb 
      par_2_max_outb = par_2_max_1sigma_outb 
      par_3_max_outb = par_3_max_1sigma_outb 
      par_4_max_outb = par_4_max_1sigma_outb 
    }
    //med
    if (pcal_fiducial_cut_lvl==1){
      par_0_min_inb = par_0_min_15sigma_inb 
      par_1_min_inb = par_1_min_15sigma_inb 
      par_2_min_inb = par_2_min_15sigma_inb
      par_3_min_inb = par_3_min_15sigma_inb 
      par_4_min_inb = par_4_min_15sigma_inb 
      par_0_max_inb = par_0_max_15sigma_inb 
      par_1_max_inb = par_1_max_15sigma_inb 
      par_2_max_inb = par_2_max_15sigma_inb 
      par_3_max_inb = par_3_max_15sigma_inb 
      par_4_max_inb = par_4_max_15sigma_inb 

      par_0_min_outb = par_0_min_15sigma_outb 
      par_1_min_outb = par_1_min_15sigma_outb 
      par_2_min_outb = par_2_min_15sigma_outb
      par_3_min_outb = par_3_min_15sigma_outb 
      par_4_min_outb = par_4_min_15sigma_outb 
      par_0_max_outb = par_0_max_15sigma_outb 
      par_1_max_outb = par_1_max_15sigma_outb 
      par_2_max_outb = par_2_max_15sigma_outb 
      par_3_max_outb = par_3_max_15sigma_outb 
      par_4_max_outb = par_4_max_15sigma_outb 
    }
    //tight
    if (pcal_fiducial_cut_lvl==1){
      par_0_min_inb = par_0_min_2sigma_inb 
      par_1_min_inb = par_1_min_2sigma_inb 
      par_2_min_inb = par_2_min_2sigma_inb
      par_3_min_inb = par_3_min_2sigma_inb 
      par_4_min_inb = par_4_min_2sigma_inb 
      par_0_max_inb = par_0_max_2sigma_inb 
      par_1_max_inb = par_1_max_2sigma_inb 
      par_2_max_inb = par_2_max_2sigma_inb 
      par_3_max_inb = par_3_max_2sigma_inb 
      par_4_max_inb = par_4_max_2sigma_inb 

      par_0_min_outb = par_0_min_2sigma_outb 
      par_1_min_outb = par_1_min_2sigma_outb 
      par_2_min_outb = par_2_min_2sigma_outb
      par_3_min_outb = par_3_min_2sigma_outb 
      par_4_min_outb = par_4_min_2sigma_outb 
      par_0_max_outb = par_0_max_2sigma_outb 
      par_1_max_outb = par_1_max_2sigma_outb 
      par_2_max_outb = par_2_max_2sigma_outb 
      par_3_max_outb = par_3_max_2sigma_outb 
      par_4_max_outb = par_4_max_2sigma_outb 
    }
  }

  void setGammaCutParameters(magnetic_field_config){
    println('[GammaFromEvent::setGammaCutParameters] -> setting gamma cut parameters for field ' + magnetic_field_config)
    if( magnetic_field_config == "outbending" ){
      println('[GammaFromEvent::setGammaCutParameters] -> setting parameters for outbending')

      par_0_min = par_0_min_outb 
      par_1_min = par_1_min_outb 
      par_2_min = par_2_min_outb
      par_3_min = par_3_min_outb 
      par_4_min = par_4_min_outb 
      par_0_max = par_0_max_outb 
      par_1_max = par_1_max_outb 
      par_2_max = par_2_max_outb 
      par_3_max = par_3_max_outb 
      par_4_max = par_4_max_outb 
    }
    else if( magnetic_field_config == "inbending" ){
      println('[GammaFromEvent::setGammaCutParameters] -> setting parameters for inbending')

      par_0_min = par_0_min_inb 
      par_1_min = par_1_min_inb 
      par_2_min = par_2_min_inb
      par_3_min = par_3_min_inb 
      par_4_min = par_4_min_inb 
      par_0_max = par_0_max_inb 
      par_1_max = par_1_max_inb 
      par_2_max = par_2_max_inb 
      par_3_max = par_3_max_inb 
      par_4_max = par_4_max_inb 
    }
  }


  def passGammaEBPIDCut = { event, index ->
    return (event.pid[index] == ebPID)
  }


  def passGammaPCALFiducialCut = { event, index ->

    if (event.status[index]<2000) return true

    if (event.pcal_status.contains(index)){   
        def sec = event.pcal_sector.get(index)-1
        def hit = event.pcal.get(index)
            if (hit){
                def (theta_PCAL, phi_PCAL) = thetaphifromhit([hit.x, hit.y, hit.z], sec)
                def (phi_PCAL_min, phi_PCAL_max) = PCALphimaxmin(theta_PCAL, sec)
                return phi_PCAL > phi_PCAL_min && phi_PCAL < phi_PCAL_max
            } else {
                return false
            }
    }
    return false
  }

  def thetaphifromhit(hit, sec){
      def theta_PCAL = Math.toDegrees(
          Math.acos(hit.get(2)/Math.sqrt(hit.get(0)*hit.get(0) + hit.get(1)*hit.get(1) + hit.get(2)*hit.get(2))))
      def phi_PCAL_raw = Math.toDegrees(Math.atan2(hit.get(1), hit.get(0)))
      def phi_PCAL = 0;
      if(sec+1 == 1) phi_PCAL = phi_PCAL_raw;
      if(sec+1 == 2) phi_PCAL = phi_PCAL_raw - 60;
      if(sec+1 == 3) phi_PCAL = phi_PCAL_raw - 120;
      if(sec+1 == 4 && phi_PCAL_raw > 0) phi_PCAL = phi_PCAL_raw - 180;
      if(sec+1 == 4 && phi_PCAL_raw < 0) phi_PCAL = phi_PCAL_raw + 180;
      if(sec+1 == 5) phi_PCAL = phi_PCAL_raw + 120;
      if(sec+1 == 6) phi_PCAL = phi_PCAL_raw + 60;
      return [theta_PCAL, phi_PCAL]
  }

  def PCALphimaxmin(theta_PCAL, sec){

    double min_phi = par_0_min[sec] + par_1_min[sec] * Math.log(theta_PCAL) + par_2_min[sec] * theta_PCAL + par_3_min[sec] *theta_PCAL*theta_PCAL + par_4_min[sec] *theta_PCAL*theta_PCAL*theta_PCAL;
    double max_phi = par_0_max[sec] + par_1_max[sec] * Math.log(theta_PCAL) + par_2_max[sec] * theta_PCAL + par_3_max[sec] *theta_PCAL*theta_PCAL + par_4_max[sec] *theta_PCAL*theta_PCAL*theta_PCAL;

    return [min_phi, max_phi];
  }

}