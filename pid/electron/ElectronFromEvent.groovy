package pid.electron

import org.jlab.detector.base.DetectorType
import event.Event

class ElectronFromEvent {

    def ebeam = 10.6
    def ebPID = 11
    def min_nphe = 2
    def nphe_loose = -1
    def nphe_med = 0
    def nphe_tight = 1
    
    //vertex wide enough for all sectors
    /// def min_vz = -12
    /// def max_vz = 9

    def min_ecal_inner_dep = 0.06
    def ecal_inner_dep_loose = -0.01
    def ecal_inner_dep_med = 0.0
    def ecal_inner_dep_tight = 0.01

    //calorimeter fiducial cut
    //def min_u = 0
    //def max_u = 420
    //def min_v = 0
    //def max_v = 420
    //def min_w = 0
    //def max_w = 420

    
    //def s_min_u_outb = [33, 26, 34, 22, 27, 25]
    //def s_max_u_outb = [408, 408, 408, 408, 408, 408]
    
    //def s_min_v_outb = [16.0, 10.5, 17.0, 14.25, 18.0, 11.0]
    //def s_max_v_outb = [400, 400, 400, 400, 400, 400]
    
    //def s_min_w_outb = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    //def s_max_w_outb = [400, 400, 400, 400, 400, 400]

    def min_u_tight_inb = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_u_med_inb   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_u_loose_inb = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_inb = [398, 398, 398, 398, 398, 398]
    def max_u_med_inb   = [408, 408, 408, 408, 408, 408]
    def max_u_loose_inb = [420, 420, 420, 420, 420, 420] 
    
    def min_u_cuts_inb = [min_u_loose_inb, min_u_med_inb, min_u_tight_inb]
    def max_u_cuts_inb = [max_u_loose_inb, max_u_med_inb, max_u_tight_inb]

    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_inb = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_v_med_inb   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_v_loose_inb = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // the maximum of v is never reached
    def max_v_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_v_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_inb = [min_v_loose_inb, min_v_med_inb, min_v_tight_inb]
    def max_v_cuts_inb = [max_v_loose_inb, max_v_med_inb, max_v_tight_inb]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_w_tight_inb = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_w_med_inb   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_w_loose_inb = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // the maximum of w is never reached
    def max_w_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_w_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_w_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_w_cuts_inb = [min_w_loose_inb, min_w_med_inb, min_w_tight_inb]
    def max_w_cuts_inb = [max_w_loose_inb, max_w_med_inb, max_w_tight_inb]
   
    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_u_tight_out = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_u_med_out   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_u_loose_out = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_out = [398, 398, 398, 398, 398, 398]
    def max_u_med_out   = [408, 408, 408, 408, 408, 408] 
    def max_u_loose_out = [420, 420, 420, 420, 420, 420] 

    def min_u_cuts_outb = [min_u_loose_out, min_u_med_out, min_u_tight_out]
    def max_u_cuts_outb = [max_u_loose_out, max_u_med_out, max_u_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_out = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_v_med_out   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_v_loose_out = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // the maximum of v is never reached
    def max_v_tight_out = [400, 400, 400, 400, 400, 400]
    def max_v_med_out   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_out = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_outb = [min_v_loose_out, min_v_med_out, min_v_tight_out]
    def max_v_cuts_outb = [max_v_loose_out, max_v_med_out, max_v_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_w_tight_out = [19.0, 19.0, 19.0, 19.0, 19.0, 19.0]
    def min_w_med_out   = [14.0, 14.0, 14.0, 14.0, 14.0, 14.0]
    def min_w_loose_out = [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
    
    // the maximum of w is never reached
    def max_w_tight_out = [400, 400, 400, 400, 400, 400]
    def max_w_med_out   = [400, 400, 400, 400, 400, 400]
    def max_w_loose_out = [400, 400, 400, 400, 400, 400]

    def min_w_cuts_outb = [min_w_loose_out, min_w_med_out, min_w_tight_out]
    def max_w_cuts_outb = [max_w_loose_out, max_w_med_out, max_w_tight_out]

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    def sigma_range = 5.0
    def p0mean_inb = [0.105631, 0.11551, 0.112799, 0.109937, 0.116249, 0.119057]
    def p1mean_inb = [-0.153951, -0.0253273, -0.125718, 0.165414, 0.0768411, 0.0555026]
    def p2mean_inb = [0.00860091, 0.00706291, 0.00908884, 0.00499666, 0.00448701, 0.00558927]
    def p3mean_inb = [-0.000821675, -0.000711488, -0.000930922, -0.000298311, -0.000455716, -0.000657084]
    def p0sigma_inb = [0.0149613, 0.0115116, 0.00580737, 0.0106817, 0.012667, 0.00553471]
    def p1sigma_inb = [0.00700773, 0.0116193, 0.0202375, 0.0126958, 0.00892239, 0.0216206]
    
    def p0mean_outb = [0.105467, 0.115261, 0.127793, 0.113359, 0.112263, 0.113507]
    def p1mean_outb = [-0.135178, 0.135808, 0.903412, 0.598274, -0.0466815, 0.0550123]
    def p2mean_outb = [0.00996842, 0.00672508, -0.00035721, 0.00470925, 0.00588451, 0.00923385]
    def p3mean_outb = [-0.000754536, -0.000515365, -0.000108273, -0.000447278, -0.000358148, -0.00074643]
    def p0sigma_outb = [0.00683747, 0.0065199, 0.00297734, 0.00759701, 0.0093309, 0.00591988]
    def p1sigma_outb = [0.0180228, 0.0183979, 0.0250332, 0.0155001, 0.0137594, 0.0215643]

    def p0mean = null
    def p1mean = null
    def p2mean = null
    def p3mean = null
    def p0sigma = null
    def p1sigma = null

    //dcr1,2,3 fiducial cut
    //require functional form
    def sect_angle_coverage = 60
    
    

    //////////////////////////////////////////////
    // use the tracking bank to get the sector for tracks

    def dcR1_loose = -1
    def dcR1_med = 0
    def dcR1_tight = 1

    def dcR2_loose = -2
    def dcR2_med = 0
    def dcR2_tight = 2

    def dcR3_loose = -3
    def dcR3_med = 0
    def dcR3_tight = 3

    def radiusR1_inb = 29
    def radiusR1_outb = 29
       
    def heightR1_inb = [27, 22, 22, 27, 22, 22]
    def heightR1_outb = [15, 15, 15, 15, 15, 15]
    
    def radiusR2_inb = 38
    def radiusR2_outb = 39

    def heightR2_inb  = [40, 34, 34, 40, 34, 34]
    def heightR2_outb = [25, 25, 25, 25, 25, 25]

    def radiusR3_inb = 50
    def radiusR3_outb = 65
    
    def heightR3_inb = [47, 39, 39, 47, 39, 39]
    def heightR3_outb = [48, 48, 48, 48, 48, 48]

    // updated DC cuts

    def reg1_min_sec1, reg1_min_sec1_inb, reg1_min_sec1_outb
    def reg1_max_sec1, reg1_max_sec1_inb, reg1_max_sec1_outb
    def reg2_min_sec1, reg2_min_sec1_inb, reg2_min_sec1_outb
    def reg2_max_sec1, reg2_max_sec1_inb, reg2_max_sec1_outb
    def reg3_min_sec1, reg3_min_sec1_inb, reg3_min_sec1_outb
    def reg3_max_sec1, reg3_max_sec1_inb, reg3_max_sec1_outb
    def reg1_min_sec2, reg1_min_sec2_inb, reg1_min_sec2_outb
    def reg1_max_sec2, reg1_max_sec2_inb, reg1_max_sec2_outb
    def reg2_min_sec2, reg2_min_sec2_inb, reg2_min_sec2_outb
    def reg2_max_sec2, reg2_max_sec2_inb, reg2_max_sec2_outb
    def reg3_min_sec2, reg3_min_sec2_inb, reg3_min_sec2_outb
    def reg3_max_sec2, reg3_max_sec2_inb, reg3_max_sec2_outb
    def reg1_min_sec3, reg1_min_sec3_inb, reg1_min_sec3_outb
    def reg1_max_sec3, reg1_max_sec3_inb, reg1_max_sec3_outb
    def reg2_min_sec3, reg2_min_sec3_inb, reg2_min_sec3_outb
    def reg2_max_sec3, reg2_max_sec3_inb, reg2_max_sec3_outb
    def reg3_min_sec3, reg3_min_sec3_inb, reg3_min_sec3_outb
    def reg3_max_sec3, reg3_max_sec3_inb, reg3_max_sec3_outb
    def reg1_min_sec4, reg1_min_sec4_inb, reg1_min_sec4_outb
    def reg1_max_sec4, reg1_max_sec4_inb, reg1_max_sec4_outb
    def reg2_min_sec4, reg2_min_sec4_inb, reg2_min_sec4_outb
    def reg2_max_sec4, reg2_max_sec4_inb, reg2_max_sec4_outb
    def reg3_min_sec4, reg3_min_sec4_inb, reg3_min_sec4_outb
    def reg3_max_sec4, reg3_max_sec4_inb, reg3_max_sec4_outb
    def reg1_min_sec5, reg1_min_sec5_inb, reg1_min_sec5_outb
    def reg1_max_sec5, reg1_max_sec5_inb, reg1_max_sec5_outb
    def reg2_min_sec5, reg2_min_sec5_inb, reg2_min_sec5_outb
    def reg2_max_sec5, reg2_max_sec5_inb, reg2_max_sec5_outb
    def reg3_min_sec5, reg3_min_sec5_inb, reg3_min_sec5_outb
    def reg3_max_sec5, reg3_max_sec5_inb, reg3_max_sec5_outb
    def reg1_min_sec6, reg1_min_sec6_inb, reg1_min_sec6_outb
    def reg1_max_sec6, reg1_max_sec6_inb, reg1_max_sec6_outb
    def reg2_min_sec6, reg2_min_sec6_inb, reg2_min_sec6_outb
    def reg2_max_sec6, reg2_max_sec6_inb, reg2_max_sec6_outb
    def reg3_min_sec6, reg3_min_sec6_inb, reg3_min_sec6_outb
    def reg3_max_sec6, reg3_max_sec6_inb, reg3_max_sec6_outb
    // inbending cut parameters (adjusted):
    // medium (50 %):

    def reg1_min_sec1_inb_med = [25.6464, -26.0938, 1.36205, 0.0083525, -0.000433074];
    def reg1_max_sec1_inb_med = [-21.002, 19.3141, 0.183817, -0.061424, 0.00117881];
    def reg2_min_sec1_inb_med = [20.8156, -19.1708, -0.244432, 0.0652386, -0.00120129];
    def reg2_max_sec1_inb_med = [-20.2319, 19.2371, 0.0491417, -0.0527963, 0.000995253];
    def reg3_min_sec1_inb_med = [22.9342, -24.6621, 1.25766, 0.0133507, -0.00049441];
    def reg3_max_sec1_inb_med = [-8.72949, 4.50811, 3.13299, -0.15452, 0.00234593];
    def reg1_min_sec2_inb_med = [16.0752, -15.8399, -0.65495, 0.0740595, -0.00133068];
    def reg1_max_sec2_inb_med = [-40.921, 40.2558, -3.67186, 0.0498104, -0.000151358];
    def reg2_min_sec2_inb_med = [20.5524, -23.7002, 1.20081, 0.0126066, -0.000501927];
    def reg2_max_sec2_inb_med = [-22.7933, 18.397, 0.707315, -0.0848617, 0.00150965];
    def reg3_min_sec2_inb_med = [7.84857, -9.62769, -1.39224, 0.0876156, -0.00137369];
    def reg3_max_sec2_inb_med = [-24.5099, 24.8373, -1.22162, -0.0135951, 0.00051455];
    def reg1_min_sec3_inb_med = [18.2067, -12.322, -1.79675, 0.111661, -0.0017831];
    def reg1_max_sec3_inb_med = [-20.4639, 11.3285, 2.56498, -0.149755, 0.00243423];
    def reg2_min_sec3_inb_med = [17.5093, -13.0897, -1.44801, 0.0995944, -0.00163973];
    def reg2_max_sec3_inb_med = [-29.4241, 25.1723, -0.643032, -0.0433457, 0.0010157];
    def reg3_min_sec3_inb_med = [19.6412, -18.3363, -0.139482, 0.0623644, -0.00122319];
    def reg3_max_sec3_inb_med = [-23.4069, 20.1087, 0.0111006, -0.0595069, 0.00123094];
    def reg1_min_sec4_inb_med = [27.7329, -24.8668, 0.855257, 0.0229504, -0.000588103];
    def reg1_max_sec4_inb_med = [-25.1595, 24.3626, -0.826633, -0.0280833, 0.000694613];
    def reg2_min_sec4_inb_med = [35.3061, -35.9815, 3.34709, -0.0587318, 0.000517383];
    def reg2_max_sec4_inb_med = [-23.5005, 23.5232, -0.872864, -0.0214562, 0.000570724];
    def reg3_min_sec4_inb_med = [13.1517, -8.73139, -2.18634, 0.115133, -0.00170354];
    def reg3_max_sec4_inb_med = [-14.1635, 13.6887, 0.949846, -0.0787818, 0.00128683];
    def reg1_min_sec5_inb_med = [23.5653, -17.285, -1.21053, 0.108096, -0.00189499];
    def reg1_max_sec5_inb_med = [-23.3017, 16.1351, 1.45328, -0.112999, 0.00193088];
    def reg2_min_sec5_inb_med = [30.5631, -30.3887, 2.11072, -0.00879513, -0.000281135];
    def reg2_max_sec5_inb_med = [-18.7762, 11.5513, 2.23634, -0.132, 0.00208225];
    def reg3_min_sec5_inb_med = [14.2375, -10.0129, -2.03418, 0.1236, -0.00199283];
    def reg3_max_sec5_inb_med = [-23.8159, 21.7783, -0.373554, -0.0463973, 0.0010055];
    def reg1_min_sec6_inb_med = [34.9766, -38.2253, 3.80455, -0.0624294, 0.000385789];
    def reg1_max_sec6_inb_med = [-23.6287, 21.0793, 0.0916929, -0.0641014, 0.00125089];
    def reg2_min_sec6_inb_med = [21.6462, -23.9106, 1.1985, 0.0108284, -0.000435214];
    def reg2_max_sec6_inb_med = [-22.1881, 21.3921, -0.325758, -0.0406509, 0.000821165];
    def reg3_min_sec6_inb_med = [12.3393, -12.5734, -1.07645, 0.0836258, -0.00137295];
    def reg3_max_sec6_inb_med = [-20.6041, 23.7403, -1.43195, 0.00507544, 0.000160285];

    // loose (35 %):
    def reg1_min_sec1_inb_loose = [24.9198, -26.6682, 1.22328, 0.0326184, 0.0326184];
    def reg1_max_sec1_inb_loose = [-38.0903, 44.8413, -5.13684, 0.0913599, -0.000566229];
    def reg2_min_sec1_inb_loose = [34.3872, -41.6892, 4.78076, -0.0904016, 0.000691157];
    def reg2_max_sec1_inb_loose = [-36.582, 45.4679, -5.70527, 0.121105, -0.00107292];
    def reg3_min_sec1_inb_loose = [38.5882, -51.7707, 7.65247, -0.198183, 0.00222671];
    def reg3_max_sec1_inb_loose = [-35.1098, 46.9821, -6.62713, 0.16512, -0.00177347];
    def reg1_min_sec2_inb_loose = [26.3964, -32.5314, 2.87032, -0.0261837, -0.0261837];
    def reg1_max_sec2_inb_loose = [-45.4528, 53.4467, -6.83895, 0.144609, -0.00127358];
    def reg2_min_sec2_inb_loose = [29.5885, -38.6103, 4.45494, -0.084114, 0.000626459];
    def reg2_max_sec2_inb_loose = [-42.6042, 50.783, -6.44381, 0.135305, -0.00117332];
    def reg3_min_sec2_inb_loose = [21.3286, -33.2316, 4.16718, -0.0947928, 0.000950229];
    def reg3_max_sec2_inb_loose = [-26.2955, 32.832, -3.21406, 0.0443669, -0.000110018];
    def reg1_min_sec3_inb_loose = [38.3683, -39.555, 3.48786, -0.0311045, -0.0311045];
    def reg1_max_sec3_inb_loose = [-44.7002, 48.5238, -5.33933, 0.0857334, -0.000389759];
    def reg2_min_sec3_inb_loose = [26.2411, -26.4956, 1.17454, 0.0314168, -0.000923764];
    def reg2_max_sec3_inb_loose = [-49.1042, 53.7936, -6.50797, 0.125435, -0.000940508];
    def reg3_min_sec3_inb_loose = [30.2087, -36.2303, 3.94272, -0.070659, 0.000484111];
    def reg3_max_sec3_inb_loose = [-25.1075, 27.5492, -1.81589, -0.00372993, 0.000534551];
    def reg1_min_sec4_inb_loose = [45.5355, -50.2561, 5.91736, -0.116203, -0.116203];
    def reg1_max_sec4_inb_loose = [-45.3956, 52.8208, -6.57583, 0.13245, -0.00106037];
    def reg2_min_sec4_inb_loose = [42.9007, -48.8609, 6.02101, -0.13073, 0.00125121];
    def reg2_max_sec4_inb_loose = [-32.1523, 38.6933, -4.15879, 0.0722104, -0.000474606];
    def reg3_min_sec4_inb_loose = [34.0991, -40.7811, 4.8997, -0.107033, 0.00105086];
    def reg3_max_sec4_inb_loose = [-39.189, 52.9641, -8.01756, 0.214872, -0.00250226];
    def reg1_min_sec5_inb_loose = [52.0273, -61.3244, 8.45187, -0.193011, -0.193011];
    def reg1_max_sec5_inb_loose = [-50.6284, 54.3673, -6.38031, 0.116905, -0.000778886];
    def reg2_min_sec5_inb_loose = [32.1495, -35.3084, 3.11177, -0.0314703, -8.39833e-05];
    def reg2_max_sec5_inb_loose = [-47.0576, 50.676, -5.80697, 0.10389, -0.000677858];
    def reg3_min_sec5_inb_loose = [27.5659, -34.8935, 3.81989, -0.0655272, 0.000389246];
    def reg3_max_sec5_inb_loose = [-25.7138, 27.0011, -1.53335, -0.0152295, 0.000694336];
    def reg1_min_sec6_inb_loose = [35.6586, -43.3351, 5.05374, -0.0961372, -0.0961372];
    def reg1_max_sec6_inb_loose = [-43.6986, 52.1413, -6.61602, 0.137247, -0.00115149];
    def reg2_min_sec6_inb_loose = [33.3586, -42.3939, 5.17582, -0.106612, 0.000910642];
    def reg2_max_sec6_inb_loose = [-38.3305, 46.3099, -5.57615, 0.109409, -0.000869664];
    def reg3_min_sec6_inb_loose = [20.9702, -29.641, 3.08715, -0.0542694, 0.00036852];
    def reg3_max_sec6_inb_loose = [-36.3081, 47.6882, -6.51571, 0.152934, -0.00153014];

    // tight (65 %):
    def reg1_min_sec1_inb_tight = [10.8127, 0.867417, -5.23315, 0.236924, 0.236924];
    def reg1_max_sec1_inb_tight = [-1.05905, -11.11, 6.99604, -0.280718, 0.00394115];
    def reg2_min_sec1_inb_tight = [15.6429, -8.51554, -2.8828, 0.158688, -0.00249913];
    def reg2_max_sec1_inb_tight = [4.94591, -19.5039, 8.81259, -0.33865, 0.00465444];
    def reg3_min_sec1_inb_tight = [1.34092, 9.28247, -6.60044, 0.281318, -0.00409396];
    def reg3_max_sec1_inb_tight = [7.95213, -22.3926, 9.38108, -0.362289, 0.00504018];
    def reg1_min_sec2_inb_tight = [6.83952, 1.15885, -4.75388, 0.212145, 0.212145];
    def reg1_max_sec2_inb_tight = [-13.5338, 2.90293, 4.04606, -0.18406, 0.00270977];
    def reg2_min_sec2_inb_tight = [2.23926, 6.12481, -5.72273, 0.242722, -0.00348261];
    def reg2_max_sec2_inb_tight = [-19.2228, 12.2586, 1.80921, -0.109461, 0.00172974];
    def reg3_min_sec2_inb_tight = [-7.44705, 15.9478, -7.51342, 0.2998, -0.00422782];
    def reg3_max_sec2_inb_tight = [-8.05779, -1.88087, 4.79989, -0.208209, 0.00301469];
    def reg1_min_sec3_inb_tight = [5.36591, 8.06724, -6.39578, 0.262518, 0.262518];
    def reg1_max_sec3_inb_tight = [-8.20889, -5.7771, 5.98507, -0.245252, 0.00348177];
    def reg2_min_sec3_inb_tight = [-2.61308, 17.4864, -8.28425, 0.321424, -0.00449555];
    def reg2_max_sec3_inb_tight = [-15.7028, 5.06188, 3.58882, -0.169598, 0.00252206];
    def reg3_min_sec3_inb_tight = [-1.47028, 12.199, -6.67515, 0.261044, -0.00361062];
    def reg3_max_sec3_inb_tight = [-2.49195, -10.9866, 6.77373, -0.267217, 0.0037212];
    def reg1_min_sec4_inb_tight = [7.28085, 4.43634, -5.4954, 0.22469, 0.22469];
    def reg1_max_sec4_inb_tight = [-2.88726, -7.75256, 6.11348, -0.246024, 0.00342];
    def reg2_min_sec4_inb_tight = [11.1628, -0.875717, -4.40181, 0.191682, -0.00269305];
    def reg2_max_sec4_inb_tight = [4.98009, -20.3121, 9.08305, -0.347362, 0.00476636];
    def reg3_min_sec4_inb_tight = [-1.53387, 16.9129, -8.4497, 0.334303, -0.00465195];
    def reg3_max_sec4_inb_tight = [5.76932, -17.8998, 8.19437, -0.317309, 0.00436422];
    def reg1_min_sec5_inb_tight = [16.2619, -7.2257, -3.02427, 0.151797, 0.151797];
    def reg1_max_sec5_inb_tight = [-8.7963, -3.03534, 5.17438, -0.214586, 0.0030239];
    def reg2_min_sec5_inb_tight = [14.656, -7.43444, -2.74998, 0.141668, -0.00216617];
    def reg2_max_sec5_inb_tight = [2.24964, -18.1672, 8.48444, -0.321566, 0.00438927];
    def reg3_min_sec5_inb_tight = [3.38717, 7.1616, -5.73249, 0.231652, -0.00320412];
    def reg3_max_sec5_inb_tight = [4.63583, -20.0558, 8.78226, -0.33313, 0.0045258];
    def reg1_min_sec6_inb_tight = [8.74274, 0.225052, -4.62671, 0.205157, 0.205157];
    def reg1_max_sec6_inb_tight = [-8.82084, 1.73358, 3.85433, -0.168173, 0.0024055];
    def reg2_min_sec6_inb_tight = [3.86614, 7.51372, -6.36287, 0.268641, -0.00385432];
    def reg2_max_sec6_inb_tight = [-7.0834, -0.809866, 4.46815, -0.18982, 0.00266762];
    def reg3_min_sec6_inb_tight = [3.12244, 3.84008, -5.05955, 0.223398, -0.00326745];
    def reg3_max_sec6_inb_tight = [5.70817, -16.9736, 7.8074, -0.295122, 0.00398425];

    // outbending cut parameters (not adjusted, same as inbending):

    //_medium (50 %):
    def reg1_min_sec1_outb_med = [25.6464, -26.0938, 1.36205, 0.0083525, 0.0083525];
    def reg1_max_sec1_outb_med = [-21.002, 19.3141, 0.183817, -0.061424, 0.00117881];
    def reg2_min_sec1_outb_med = [20.8156, -19.1708, -0.244432, 0.0652386, -0.00120129];
    def reg2_max_sec1_outb_med = [-20.2319, 19.2371, 0.0491417, -0.0527963, 0.000995253];
    def reg3_min_sec1_outb_med = [22.9342, -24.6621, 1.25766, 0.0133507, -0.00049441];
    def reg3_max_sec1_outb_med = [-8.72949, 4.50811, 3.13299, -0.15452, 0.00234593];
    def reg1_min_sec2_outb_med = [16.0752, -15.8399, -0.65495, 0.0740595, 0.0740595];
    def reg1_max_sec2_outb_med = [-40.921, 40.2558, -3.67186, 0.0498104, -0.000151358];
    def reg2_min_sec2_outb_med = [20.5524, -23.7002, 1.20081, 0.0126066, -0.000501927];
    def reg2_max_sec2_outb_med = [-22.7933, 18.397, 0.707315, -0.0848617, 0.00150965];
    def reg3_min_sec2_outb_med = [7.84857, -9.62769, -1.39224, 0.0876156, -0.00137369];
    def reg3_max_sec2_outb_med = [-24.5099, 24.8373, -1.22162, -0.0135951, 0.00051455];
    def reg1_min_sec3_outb_med = [18.2067, -12.322, -1.79675, 0.111661, 0.111661];
    def reg1_max_sec3_outb_med = [-20.4639, 11.3285, 2.56498, -0.149755, 0.00243423];
    def reg2_min_sec3_outb_med = [17.5093, -13.0897, -1.44801, 0.0995944, -0.00163973];
    def reg2_max_sec3_outb_med = [-29.4241, 25.1723, -0.643032, -0.0433457, 0.0010157];
    def reg3_min_sec3_outb_med = [19.6412, -18.3363, -0.139482, 0.0623644, -0.00122319];
    def reg3_max_sec3_outb_med = [-23.4069, 20.1087, 0.0111006, -0.0595069, 0.00123094];
    def reg1_min_sec4_outb_med = [27.7329, -24.8668, 0.855257, 0.0229504, 0.0229504];
    def reg1_max_sec4_outb_med = [-25.1595, 24.3626, -0.826633, -0.0280833, 0.000694613];
    def reg2_min_sec4_outb_med = [35.3061, -35.9815, 3.34709, -0.0587318, 0.000517383];
    def reg2_max_sec4_outb_med = [-23.5005, 23.5232, -0.872864, -0.0214562, 0.000570724];
    def reg3_min_sec4_outb_med = [13.1517, -8.73139, -2.18634, 0.115133, -0.00170354];
    def reg3_max_sec4_outb_med = [-14.1635, 13.6887, 0.949846, -0.0787818, 0.00128683];
    def reg1_min_sec5_outb_med = [23.5653, -17.285, -1.21053, 0.108096, 0.108096];
    def reg1_max_sec5_outb_med = [-23.3017, 16.1351, 1.45328, -0.112999, 0.00193088];
    def reg2_min_sec5_outb_med = [30.5631, -30.3887, 2.11072, -0.00879513, -0.000281135];
    def reg2_max_sec5_outb_med = [-18.7762, 11.5513, 2.23634, -0.132, 0.00208225];
    def reg3_min_sec5_outb_med = [14.2375, -10.0129, -2.03418, 0.1236, -0.00199283];
    def reg3_max_sec5_outb_med = [-23.8159, 21.7783, -0.373554, -0.0463973, 0.0010055];
    def reg1_min_sec6_outb_med = [34.9766, -38.2253, 3.80455, -0.0624294, -0.0624294];
    def reg1_max_sec6_outb_med = [-23.6287, 21.0793, 0.0916929, -0.0641014, 0.00125089];
    def reg2_min_sec6_outb_med = [21.6462, -23.9106, 1.1985, 0.0108284, -0.000435214];
    def reg2_max_sec6_outb_med = [-22.1881, 21.3921, -0.325758, -0.0406509, 0.000821165];
    def reg3_min_sec6_outb_med = [12.3393, -12.5734, -1.07645, 0.0836258, -0.00137295];
    def reg3_max_sec6_outb_med = [-20.6041, 23.7403, -1.43195, 0.00507544, 0.000160285];

    // loose (35 %):
    def reg1_min_sec1_outb_loose = [24.9198, -26.6682, 1.22328, 0.0326184, 0.0326184];
    def reg1_max_sec1_outb_loose = [-38.0903, 44.8413, -5.13684, 0.0913599, -0.000566229];
    def reg2_min_sec1_outb_loose = [34.3872, -41.6892, 4.78076, -0.0904016, 0.000691157];
    def reg2_max_sec1_outb_loose = [-36.582, 45.4679, -5.70527, 0.121105, -0.00107292];
    def reg3_min_sec1_outb_loose = [38.5882, -51.7707, 7.65247, -0.198183, 0.00222671];
    def reg3_max_sec1_outb_loose = [-35.1098, 46.9821, -6.62713, 0.16512, -0.00177347];
    def reg1_min_sec2_outb_loose = [26.3964, -32.5314, 2.87032, -0.0261837, -0.0261837];
    def reg1_max_sec2_outb_loose = [-45.4528, 53.4467, -6.83895, 0.144609, -0.00127358];
    def reg2_min_sec2_outb_loose = [29.5885, -38.6103, 4.45494, -0.084114, 0.000626459];
    def reg2_max_sec2_outb_loose = [-42.6042, 50.783, -6.44381, 0.135305, -0.00117332];
    def reg3_min_sec2_outb_loose = [21.3286, -33.2316, 4.16718, -0.0947928, 0.000950229];
    def reg3_max_sec2_outb_loose = [-26.2955, 32.832, -3.21406, 0.0443669, -0.000110018];
    def reg1_min_sec3_outb_loose = [38.3683, -39.555, 3.48786, -0.0311045, -0.0311045];
    def reg1_max_sec3_outb_loose = [-44.7002, 48.5238, -5.33933, 0.0857334, -0.000389759];
    def reg2_min_sec3_outb_loose = [26.2411, -26.4956, 1.17454, 0.0314168, -0.000923764];
    def reg2_max_sec3_outb_loose = [-49.1042, 53.7936, -6.50797, 0.125435, -0.000940508];
    def reg3_min_sec3_outb_loose = [30.2087, -36.2303, 3.94272, -0.070659, 0.000484111];
    def reg3_max_sec3_outb_loose = [-25.1075, 27.5492, -1.81589, -0.00372993, 0.000534551];
    def reg1_min_sec4_outb_loose = [45.5355, -50.2561, 5.91736, -0.116203, -0.116203];
    def reg1_max_sec4_outb_loose = [-45.3956, 52.8208, -6.57583, 0.13245, -0.00106037];
    def reg2_min_sec4_outb_loose = [42.9007, -48.8609, 6.02101, -0.13073, 0.00125121];
    def reg2_max_sec4_outb_loose = [-32.1523, 38.6933, -4.15879, 0.0722104, -0.000474606];
    def reg3_min_sec4_outb_loose = [34.0991, -40.7811, 4.8997, -0.107033, 0.00105086];
    def reg3_max_sec4_outb_loose = [-39.189, 52.9641, -8.01756, 0.214872, -0.00250226];
    def reg1_min_sec5_outb_loose = [52.0273, -61.3244, 8.45187, -0.193011, -0.193011];
    def reg1_max_sec5_outb_loose = [-50.6284, 54.3673, -6.38031, 0.116905, -0.000778886];
    def reg2_min_sec5_outb_loose = [32.1495, -35.3084, 3.11177, -0.0314703, -8.39833e-05];
    def reg2_max_sec5_outb_loose = [-47.0576, 50.676, -5.80697, 0.10389, -0.000677858];
    def reg3_min_sec5_outb_loose = [27.5659, -34.8935, 3.81989, -0.0655272, 0.000389246];
    def reg3_max_sec5_outb_loose = [-25.7138, 27.0011, -1.53335, -0.0152295, 0.000694336];
    def reg1_min_sec6_outb_loose = [35.6586, -43.3351, 5.05374, -0.0961372, -0.0961372];
    def reg1_max_sec6_outb_loose = [-43.6986, 52.1413, -6.61602, 0.137247, -0.00115149];
    def reg2_min_sec6_outb_loose = [33.3586, -42.3939, 5.17582, -0.106612, 0.000910642];
    def reg2_max_sec6_outb_loose = [-38.3305, 46.3099, -5.57615, 0.109409, -0.000869664];
    def reg3_min_sec6_outb_loose = [20.9702, -29.641, 3.08715, -0.0542694, 0.00036852];
    def reg3_max_sec6_outb_loose = [-36.3081, 47.6882, -6.51571, 0.152934, -0.00153014];

    // tight (65 %):
    def reg1_min_sec1_outb_tight = [10.8127, 0.867417, -5.23315, 0.236924, 0.236924];
    def reg1_max_sec1_outb_tight = [-1.05905, -11.11, 6.99604, -0.280718, 0.00394115];
    def reg2_min_sec1_outb_tight = [15.6429, -8.51554, -2.8828, 0.158688, -0.00249913];
    def reg2_max_sec1_outb_tight = [4.94591, -19.5039, 8.81259, -0.33865, 0.00465444];
    def reg3_min_sec1_outb_tight = [1.34092, 9.28247, -6.60044, 0.281318, -0.00409396];
    def reg3_max_sec1_outb_tight = [7.95213, -22.3926, 9.38108, -0.362289, 0.00504018];
    def reg1_min_sec2_outb_tight = [6.83952, 1.15885, -4.75388, 0.212145, 0.212145];
    def reg1_max_sec2_outb_tight = [-13.5338, 2.90293, 4.04606, -0.18406, 0.00270977];
    def reg2_min_sec2_outb_tight = [2.23926, 6.12481, -5.72273, 0.242722, -0.00348261];
    def reg2_max_sec2_outb_tight = [-19.2228, 12.2586, 1.80921, -0.109461, 0.00172974];
    def reg3_min_sec2_outb_tight = [-7.44705, 15.9478, -7.51342, 0.2998, -0.00422782];
    def reg3_max_sec2_outb_tight = [-8.05779, -1.88087, 4.79989, -0.208209, 0.00301469];
    def reg1_min_sec3_outb_tight = [5.36591, 8.06724, -6.39578, 0.262518, 0.262518];
    def reg1_max_sec3_outb_tight = [-8.20889, -5.7771, 5.98507, -0.245252, 0.00348177];
    def reg2_min_sec3_outb_tight = [-2.61308, 17.4864, -8.28425, 0.321424, -0.00449555];
    def reg2_max_sec3_outb_tight = [-15.7028, 5.06188, 3.58882, -0.169598, 0.00252206];
    def reg3_min_sec3_outb_tight = [-1.47028, 12.199, -6.67515, 0.261044, -0.00361062];
    def reg3_max_sec3_outb_tight = [-2.49195, -10.9866, 6.77373, -0.267217, 0.0037212];
    def reg1_min_sec4_outb_tight = [7.28085, 4.43634, -5.4954, 0.22469, 0.22469];
    def reg1_max_sec4_outb_tight = [-2.88726, -7.75256, 6.11348, -0.246024, 0.00342];
    def reg2_min_sec4_outb_tight = [11.1628, -0.875717, -4.40181, 0.191682, -0.00269305];
    def reg2_max_sec4_outb_tight = [4.98009, -20.3121, 9.08305, -0.347362, 0.00476636];
    def reg3_min_sec4_outb_tight = [-1.53387, 16.9129, -8.4497, 0.334303, -0.00465195];
    def reg3_max_sec4_outb_tight = [5.76932, -17.8998, 8.19437, -0.317309, 0.00436422];
    def reg1_min_sec5_outb_tight = [16.2619, -7.2257, -3.02427, 0.151797, 0.151797];
    def reg1_max_sec5_outb_tight = [-8.7963, -3.03534, 5.17438, -0.214586, 0.0030239];
    def reg2_min_sec5_outb_tight = [14.656, -7.43444, -2.74998, 0.141668, -0.00216617];
    def reg2_max_sec5_outb_tight = [2.24964, -18.1672, 8.48444, -0.321566, 0.00438927];
    def reg3_min_sec5_outb_tight = [3.38717, 7.1616, -5.73249, 0.231652, -0.00320412];
    def reg3_max_sec5_outb_tight = [4.63583, -20.0558, 8.78226, -0.33313, 0.0045258];
    def reg1_min_sec6_outb_tight = [8.74274, 0.225052, -4.62671, 0.205157, 0.205157];
    def reg1_max_sec6_outb_tight = [-8.82084, 1.73358, 3.85433, -0.168173, 0.0024055];
    def reg2_min_sec6_outb_tight = [3.86614, 7.51372, -6.36287, 0.268641, -0.00385432];
    def reg2_max_sec6_outb_tight = [-7.0834, -0.809866, 4.46815, -0.18982, 0.00266762];
    def reg3_min_sec6_outb_tight = [3.12244, 3.84008, -5.05955, 0.223398, -0.00326745];
    def reg3_max_sec6_outb_tight = [5.70817, -16.9736, 7.8074, -0.295122, 0.00398425];

    def p0_min = 0; def p1_min = 0; def p2_min = 0; def p3_min = 0; def p4_min = 0;
    def p0_max = 0; def p1_max = 0; def p2_max = 0; def p3_max = 0; def p4_max = 0;


    // vertex position cuts
    def vz_tight = -1 
    def vz_med  = 0
    def vz_loose = 1

    def vz_min_sect_inb = [-12, -12, -12, -12, -12, -12]
    def vz_max_sect_inb = [9, 9, 9, 9, 9, 9]
    
    def vz_min_sect_outb = [-14, -14, -14, -14, -14, -14]
    def vz_max_sect_outb = [5, 5, 5, 5, 5, 5]

    def s_min_u = null
    def s_max_u = null
    def s_min_v = null
    def s_max_v = null
    def s_min_w = null
    def s_max_w = null
      
    def heightR1=null
    def heightR2=null
    def heightR3=null

    def radiusR1=null
    def radiusR2=null
    def radiusR3=null

    def min_vz=null
    def max_vz=null
    def p_min=null

    def el_cut_strictness_lvl=null
    def ecal_cut_lvl = null
    def nphe_cut_lvl = null
    def vz_cut_lvl = null
    def min_u_cut_lvl = null
    def min_v_cut_lvl = null
    def min_w_cut_lvl = null
    def max_u_cut_lvl = null
    def max_v_cut_lvl = null
    def max_w_cut_lvl = null
    def dcr1_cut_lvl = null
    def dcr2_cut_lvl = null
    def dcr3_cut_lvl = null

    def anti_pion_threshold=null
            
    void setElectronCutStrictness(el_cut_strictness){
	el_cut_strictness_lvl=el_cut_strictness
	ecal_cut_lvl=el_cut_strictness["ecal_cut_lvl"]
	nphe_cut_lvl=el_cut_strictness["nphe_cut_lvl"]
	vz_cut_lvl=el_cut_strictness["vz_cut_lvl"]
	min_u_cut_lvl=el_cut_strictness["min_u_cut_lvl"]
	min_v_cut_lvl=el_cut_strictness["min_v_cut_lvl"]
	min_w_cut_lvl=el_cut_strictness["min_w_cut_lvl"]
	max_u_cut_lvl=el_cut_strictness["max_u_cut_lvl"]
	max_v_cut_lvl=el_cut_strictness["max_v_cut_lvl"]
	max_w_cut_lvl=el_cut_strictness["max_w_cut_lvl"]
	dcr1_cut_lvl=el_cut_strictness["dcr1_cut_lvl"]
	dcr2_cut_lvl=el_cut_strictness["dcr2_cut_lvl"]
	dcr3_cut_lvl=el_cut_strictness["dcr3_cut_lvl"]
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> el_cut_strictness " + el_cut_strictness)
	println("[ElectronFromEvent::setElectronCutStrictness] -> ecal_cut_lvl " + ecal_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> nphe_cut_lvl " + nphe_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> vz_cut_lvl " + vz_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_u_cut_lvl " + min_u_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_v_cut_lvl " + min_v_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> min_w_cut_lvl " + min_w_cut_lvl)
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_u_cut_lvl " + max_u_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_v_cut_lvl " + max_v_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> max_w_cut_lvl " + max_w_cut_lvl)
	
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr1_cut_lvl " + dcr1_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr2_cut_lvl " + dcr2_cut_lvl)
	println("[ElectronFromEvent::setElectronCutStrictness] -> dcr3_cut_lvl " + dcr3_cut_lvl)
	
	if( ecal_cut_lvl == 0 ) { min_ecal_inner_dep += ecal_inner_dep_loose }
	else if( ecal_cut_lvl == 1 ) { min_ecal_inner_dep += ecal_inner_dep_med }
	else if( ecal_cut_lvl == 2 ) { min_ecal_inner_dep += ecal_inner_dep_tight }

	if( nphe_cut_lvl == 0 ){ min_nphe += nphe_loose }
	else if( nphe_cut_lvl == 1 ){ min_nphe += nphe_med }
	else if( nphe_cut_lvl == 2 ){ min_nphe += nphe_tight }

	if( vz_cut_lvl == 0 ) { 
	    vz_min_sect_inb = vz_min_sect_inb.collect{ it - vz_loose }
	    vz_max_sect_inb = vz_max_sect_inb.collect{ it + vz_loose } 

	    vz_min_sect_outb = vz_min_sect_outb.collect{ it - vz_loose }
	    vz_max_sect_outb = vz_max_sect_outb.collect{ it + vz_loose } 
	}
	else if( vz_cut_lvl == 1 ) { 
	    vz_min_sect_inb = vz_min_sect_inb.collect{ it - vz_med }
	    vz_max_sect_inb = vz_max_sect_inb.collect{ it + vz_med } 

	    vz_min_sect_outb = vz_min_sect_outb.collect{ it - vz_med }
	    vz_max_sect_outb = vz_max_sect_outb.collect{ it + vz_med }  
	}
	else if( vz_cut_lvl == 2 ) { 
	    vz_min_sect_inb = vz_min_sect_inb.collect{ it - vz_tight }
	    vz_max_sect_inb = vz_max_sect_inb.collect{ it + vz_tight } 

	    vz_min_sect_outb = vz_min_sect_outb.collect{ it - vz_tight }
	    vz_max_sect_outb = vz_max_sect_outb.collect{ it + vz_tight } 
	}

	if( dcr1_cut_lvl == 0 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_loose }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_loose }
	    radiusR1_inb+=dcR1_loose
	    radiusR1_outb+=dcR1_loose

        reg1_min_sec1_inb = reg1_min_sec1_inb_loose
        reg1_max_sec1_inb = reg1_max_sec1_inb_loose
        reg1_min_sec2_inb = reg1_min_sec2_inb_loose
        reg1_max_sec2_inb = reg1_max_sec2_inb_loose
        reg1_min_sec3_inb = reg1_min_sec3_inb_loose
        reg1_max_sec3_inb = reg1_max_sec3_inb_loose
        reg1_min_sec4_inb = reg1_min_sec4_inb_loose
        reg1_max_sec4_inb = reg1_max_sec4_inb_loose
        reg1_min_sec5_inb = reg1_min_sec5_inb_loose
        reg1_max_sec5_inb = reg1_max_sec5_inb_loose
        reg1_min_sec6_inb = reg1_min_sec6_inb_loose
        reg1_max_sec6_inb = reg1_max_sec6_inb_loose

        reg1_min_sec1_outb = reg1_min_sec1_outb_loose
        reg1_max_sec1_outb = reg1_max_sec1_outb_loose
        reg1_min_sec2_outb = reg1_min_sec2_outb_loose
        reg1_max_sec2_outb = reg1_max_sec2_outb_loose
        reg1_min_sec3_outb = reg1_min_sec3_outb_loose
        reg1_max_sec3_outb = reg1_max_sec3_outb_loose
        reg1_min_sec4_outb = reg1_min_sec4_outb_loose
        reg1_max_sec4_outb = reg1_max_sec4_outb_loose
        reg1_min_sec5_outb = reg1_min_sec5_outb_loose
        reg1_max_sec5_outb = reg1_max_sec5_outb_loose
        reg1_min_sec6_outb = reg1_min_sec6_outb_loose
        reg1_max_sec6_outb = reg1_max_sec6_outb_loose

	}
	else if( dcr1_cut_lvl == 1 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_med }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_med }
	    radiusR1_inb+=dcR1_med
	    radiusR1_outb+=dcR1_med

        reg1_min_sec1_inb = reg1_min_sec1_inb_med
        reg1_max_sec1_inb = reg1_max_sec1_inb_med
        reg1_min_sec2_inb = reg1_min_sec2_inb_med
        reg1_max_sec2_inb = reg1_max_sec2_inb_med
        reg1_min_sec3_inb = reg1_min_sec3_inb_med
        reg1_max_sec3_inb = reg1_max_sec3_inb_med
        reg1_min_sec4_inb = reg1_min_sec4_inb_med
        reg1_max_sec4_inb = reg1_max_sec4_inb_med
        reg1_min_sec5_inb = reg1_min_sec5_inb_med
        reg1_max_sec5_inb = reg1_max_sec5_inb_med
        reg1_min_sec6_inb = reg1_min_sec6_inb_med
        reg1_max_sec6_inb = reg1_max_sec6_inb_med

        reg1_min_sec1_outb = reg1_min_sec1_outb_med
        reg1_max_sec1_outb = reg1_max_sec1_outb_med
        reg1_min_sec2_outb = reg1_min_sec2_outb_med
        reg1_max_sec2_outb = reg1_max_sec2_outb_med
        reg1_min_sec3_outb = reg1_min_sec3_outb_med
        reg1_max_sec3_outb = reg1_max_sec3_outb_med
        reg1_min_sec4_outb = reg1_min_sec4_outb_med
        reg1_max_sec4_outb = reg1_max_sec4_outb_med
        reg1_min_sec5_outb = reg1_min_sec5_outb_med
        reg1_max_sec5_outb = reg1_max_sec5_outb_med
        reg1_min_sec6_outb = reg1_min_sec6_outb_med
        reg1_max_sec6_outb = reg1_max_sec6_outb_med

	}
	else if( dcr1_cut_lvl == 2 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_tight }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_tight }
	    radiusR1_inb+=dcR1_tight
	    radiusR1_outb+=dcR1_tight

        reg1_min_sec1_inb = reg1_min_sec1_inb_tight
        reg1_max_sec1_inb = reg1_max_sec1_inb_tight
        reg1_min_sec2_inb = reg1_min_sec2_inb_tight
        reg1_max_sec2_inb = reg1_max_sec2_inb_tight
        reg1_min_sec3_inb = reg1_min_sec3_inb_tight
        reg1_max_sec3_inb = reg1_max_sec3_inb_tight
        reg1_min_sec4_inb = reg1_min_sec4_inb_tight
        reg1_max_sec4_inb = reg1_max_sec4_inb_tight
        reg1_min_sec5_inb = reg1_min_sec5_inb_tight
        reg1_max_sec5_inb = reg1_max_sec5_inb_tight
        reg1_min_sec6_inb = reg1_min_sec6_inb_tight
        reg1_max_sec6_inb = reg1_max_sec6_inb_tight

        reg1_min_sec1_outb = reg1_min_sec1_outb_tight
        reg1_max_sec1_outb = reg1_max_sec1_outb_tight
        reg1_min_sec2_outb = reg1_min_sec2_outb_tight
        reg1_max_sec2_outb = reg1_max_sec2_outb_tight
        reg1_min_sec3_outb = reg1_min_sec3_outb_tight
        reg1_max_sec3_outb = reg1_max_sec3_outb_tight
        reg1_min_sec4_outb = reg1_min_sec4_outb_tight
        reg1_max_sec4_outb = reg1_max_sec4_outb_tight
        reg1_min_sec5_outb = reg1_min_sec5_outb_tight
        reg1_max_sec5_outb = reg1_max_sec5_outb_tight
        reg1_min_sec6_outb = reg1_min_sec6_outb_tight
        reg1_max_sec6_outb = reg1_max_sec6_outb_tight

	}
	if( dcr2_cut_lvl == 0 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_loose }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_loose }
	    radiusR2_inb+=dcR2_loose
	    radiusR2_outb+=dcR2_loose

        reg2_min_sec1_inb = reg2_min_sec1_inb_loose
        reg2_max_sec1_inb = reg2_max_sec1_inb_loose
        reg2_min_sec2_inb = reg2_min_sec2_inb_loose
        reg2_max_sec2_inb = reg2_max_sec2_inb_loose
        reg2_min_sec3_inb = reg2_min_sec3_inb_loose
        reg2_max_sec3_inb = reg2_max_sec3_inb_loose
        reg2_min_sec4_inb = reg2_min_sec4_inb_loose
        reg2_max_sec4_inb = reg2_max_sec4_inb_loose
        reg2_min_sec5_inb = reg2_min_sec5_inb_loose
        reg2_max_sec5_inb = reg2_max_sec5_inb_loose
        reg2_min_sec6_inb = reg2_min_sec6_inb_loose
        reg2_max_sec6_inb = reg2_max_sec6_inb_loose

        reg2_min_sec1_outb = reg2_min_sec1_outb_loose
        reg2_max_sec1_outb = reg2_max_sec1_outb_loose
        reg2_min_sec2_outb = reg2_min_sec2_outb_loose
        reg2_max_sec2_outb = reg2_max_sec2_outb_loose
        reg2_min_sec3_outb = reg2_min_sec3_outb_loose
        reg2_max_sec3_outb = reg2_max_sec3_outb_loose
        reg2_min_sec4_outb = reg2_min_sec4_outb_loose
        reg2_max_sec4_outb = reg2_max_sec4_outb_loose
        reg2_min_sec5_outb = reg2_min_sec5_outb_loose
        reg2_max_sec5_outb = reg2_max_sec5_outb_loose
        reg2_min_sec6_outb = reg2_min_sec6_outb_loose
        reg2_max_sec6_outb = reg2_max_sec6_outb_loose

	}
	else if( dcr2_cut_lvl == 1 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_med }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_med }
	    radiusR2_inb+=dcR2_med
	    radiusR2_outb+=dcR2_med

        reg2_min_sec1_inb = reg2_min_sec1_inb_med
        reg2_max_sec1_inb = reg2_max_sec1_inb_med
        reg2_min_sec2_inb = reg2_min_sec2_inb_med
        reg2_max_sec2_inb = reg2_max_sec2_inb_med
        reg2_min_sec3_inb = reg2_min_sec3_inb_med
        reg2_max_sec3_inb = reg2_max_sec3_inb_med
        reg2_min_sec4_inb = reg2_min_sec4_inb_med
        reg2_max_sec4_inb = reg2_max_sec4_inb_med
        reg2_min_sec5_inb = reg2_min_sec5_inb_med
        reg2_max_sec5_inb = reg2_max_sec5_inb_med
        reg2_min_sec6_inb = reg2_min_sec6_inb_med
        reg2_max_sec6_inb = reg2_max_sec6_inb_med

        reg2_min_sec1_outb = reg2_min_sec1_outb_med
        reg2_max_sec1_outb = reg2_max_sec1_outb_med
        reg2_min_sec2_outb = reg2_min_sec2_outb_med
        reg2_max_sec2_outb = reg2_max_sec2_outb_med
        reg2_min_sec3_outb = reg2_min_sec3_outb_med
        reg2_max_sec3_outb = reg2_max_sec3_outb_med
        reg2_min_sec4_outb = reg2_min_sec4_outb_med
        reg2_max_sec4_outb = reg2_max_sec4_outb_med
        reg2_min_sec5_outb = reg2_min_sec5_outb_med
        reg2_max_sec5_outb = reg2_max_sec5_outb_med
        reg2_min_sec6_outb = reg2_min_sec6_outb_med
        reg2_max_sec6_outb = reg2_max_sec6_outb_med

	}
	else if( dcr2_cut_lvl == 2 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_tight }
	     heightR2_outb = heightR2_outb.collect{it + dcR2_tight }
	    radiusR2_inb+=dcR2_tight
	    radiusR2_outb+=dcR2_tight

        reg2_min_sec1_inb = reg2_min_sec1_inb_tight
        reg2_max_sec1_inb = reg2_max_sec1_inb_tight
        reg2_min_sec2_inb = reg2_min_sec2_inb_tight
        reg2_max_sec2_inb = reg2_max_sec2_inb_tight
        reg2_min_sec3_inb = reg2_min_sec3_inb_tight
        reg2_max_sec3_inb = reg2_max_sec3_inb_tight
        reg2_min_sec4_inb = reg2_min_sec4_inb_tight
        reg2_max_sec4_inb = reg2_max_sec4_inb_tight
        reg2_min_sec5_inb = reg2_min_sec5_inb_tight
        reg2_max_sec5_inb = reg2_max_sec5_inb_tight
        reg2_min_sec6_inb = reg2_min_sec6_inb_tight
        reg2_max_sec6_inb = reg2_max_sec6_inb_tight

        reg2_min_sec1_outb = reg2_min_sec1_outb_tight
        reg2_max_sec1_outb = reg2_max_sec1_outb_tight
        reg2_min_sec2_outb = reg2_min_sec2_outb_tight
        reg2_max_sec2_outb = reg2_max_sec2_outb_tight
        reg2_min_sec3_outb = reg2_min_sec3_outb_tight
        reg2_max_sec3_outb = reg2_max_sec3_outb_tight
        reg2_min_sec4_outb = reg2_min_sec4_outb_tight
        reg2_max_sec4_outb = reg2_max_sec4_outb_tight
        reg2_min_sec5_outb = reg2_min_sec5_outb_tight
        reg2_max_sec5_outb = reg2_max_sec5_outb_tight
        reg2_min_sec6_outb = reg2_min_sec6_outb_tight
        reg2_max_sec6_outb = reg2_max_sec6_outb_tight

	}

	if( dcr3_cut_lvl == 0 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_loose }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_loose }
	    radiusR3_inb+=dcR3_loose
	    radiusR3_outb+=dcR3_loose

        reg3_min_sec1_inb = reg3_min_sec1_inb_loose
        reg3_max_sec1_inb = reg3_max_sec1_inb_loose
        reg3_min_sec2_inb = reg3_min_sec2_inb_loose
        reg3_max_sec2_inb = reg3_max_sec2_inb_loose
        reg3_min_sec3_inb = reg3_min_sec3_inb_loose
        reg3_max_sec3_inb = reg3_max_sec3_inb_loose
        reg3_min_sec4_inb = reg3_min_sec4_inb_loose
        reg3_max_sec4_inb = reg3_max_sec4_inb_loose
        reg3_min_sec5_inb = reg3_min_sec5_inb_loose
        reg3_max_sec5_inb = reg3_max_sec5_inb_loose
        reg3_min_sec6_inb = reg3_min_sec6_inb_loose
        reg3_max_sec6_inb = reg3_max_sec6_inb_loose

        reg3_min_sec1_outb = reg3_min_sec1_outb_loose
        reg3_max_sec1_outb = reg3_max_sec1_outb_loose
        reg3_min_sec2_outb = reg3_min_sec2_outb_loose
        reg3_max_sec2_outb = reg3_max_sec2_outb_loose
        reg3_min_sec3_outb = reg3_min_sec3_outb_loose
        reg3_max_sec3_outb = reg3_max_sec3_outb_loose
        reg3_min_sec4_outb = reg3_min_sec4_outb_loose
        reg3_max_sec4_outb = reg3_max_sec4_outb_loose
        reg3_min_sec5_outb = reg3_min_sec5_outb_loose
        reg3_max_sec5_outb = reg3_max_sec5_outb_loose
        reg3_min_sec6_outb = reg3_min_sec6_outb_loose
        reg3_max_sec6_outb = reg3_max_sec6_outb_loose

	}
	else if( dcr3_cut_lvl == 1 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_med }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_med }
	    radiusR3_inb+=dcR3_med
	    radiusR3_outb+=dcR3_med

        reg3_min_sec1_inb = reg3_min_sec1_inb_med
        reg3_max_sec1_inb = reg3_max_sec1_inb_med
        reg3_min_sec2_inb = reg3_min_sec2_inb_med
        reg3_max_sec2_inb = reg3_max_sec2_inb_med
        reg3_min_sec3_inb = reg3_min_sec3_inb_med
        reg3_max_sec3_inb = reg3_max_sec3_inb_med
        reg3_min_sec4_inb = reg3_min_sec4_inb_med
        reg3_max_sec4_inb = reg3_max_sec4_inb_med
        reg3_min_sec5_inb = reg3_min_sec5_inb_med
        reg3_max_sec5_inb = reg3_max_sec5_inb_med
        reg3_min_sec6_inb = reg3_min_sec6_inb_med
        reg3_max_sec6_inb = reg3_max_sec6_inb_med

        reg3_min_sec1_outb = reg3_min_sec1_outb_med
        reg3_max_sec1_outb = reg3_max_sec1_outb_med
        reg3_min_sec2_outb = reg3_min_sec2_outb_med
        reg3_max_sec2_outb = reg3_max_sec2_outb_med
        reg3_min_sec3_outb = reg3_min_sec3_outb_med
        reg3_max_sec3_outb = reg3_max_sec3_outb_med
        reg3_min_sec4_outb = reg3_min_sec4_outb_med
        reg3_max_sec4_outb = reg3_max_sec4_outb_med
        reg3_min_sec5_outb = reg3_min_sec5_outb_med
        reg3_max_sec5_outb = reg3_max_sec5_outb_med
        reg3_min_sec6_outb = reg3_min_sec6_outb_med
        reg3_max_sec6_outb = reg3_max_sec6_outb_med

	}
	else if( dcr3_cut_lvl == 2 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_tight }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_tight }
	    radiusR3_inb+=dcR3_tight
	    radiusR3_outb+=dcR3_tight

        reg3_min_sec1_inb = reg3_min_sec1_inb_tight
        reg3_max_sec1_inb = reg3_max_sec1_inb_tight
        reg3_min_sec2_inb = reg3_min_sec2_inb_tight
        reg3_max_sec2_inb = reg3_max_sec2_inb_tight
        reg3_min_sec3_inb = reg3_min_sec3_inb_tight
        reg3_max_sec3_inb = reg3_max_sec3_inb_tight
        reg3_min_sec4_inb = reg3_min_sec4_inb_tight
        reg3_max_sec4_inb = reg3_max_sec4_inb_tight
        reg3_min_sec5_inb = reg3_min_sec5_inb_tight
        reg3_max_sec5_inb = reg3_max_sec5_inb_tight
        reg3_min_sec6_inb = reg3_min_sec6_inb_tight
        reg3_max_sec6_inb = reg3_max_sec6_inb_tight

        reg3_min_sec1_outb = reg3_min_sec1_outb_tight
        reg3_max_sec1_outb = reg3_max_sec1_outb_tight
        reg3_min_sec2_outb = reg3_min_sec2_outb_tight
        reg3_max_sec2_outb = reg3_max_sec2_outb_tight
        reg3_min_sec3_outb = reg3_min_sec3_outb_tight
        reg3_max_sec3_outb = reg3_max_sec3_outb_tight
        reg3_min_sec4_outb = reg3_min_sec4_outb_tight
        reg3_max_sec4_outb = reg3_max_sec4_outb_tight
        reg3_min_sec5_outb = reg3_min_sec5_outb_tight
        reg3_max_sec5_outb = reg3_max_sec5_outb_tight
        reg3_min_sec6_outb = reg3_min_sec6_outb_tight
        reg3_max_sec6_outb = reg3_max_sec6_outb_tight

	}	    	
	
    }

    void setElectronCutParameters(magnetic_field_config){
	println('[ElectronFromEvent::setElectronCutParameters] -> setting electron cut parameters for field ' + magnetic_field_config)
	if( magnetic_field_config == "outbending" ){
	    println('[ElectronFromEvent::setElectronCutParameters] -> setting parameters for outbending')
	    s_min_u=min_u_cuts_outb[min_u_cut_lvl] //min_u_med_out
	    s_min_v=min_v_cuts_outb[min_v_cut_lvl]
	    s_min_w=min_w_cuts_outb[min_w_cut_lvl]
	    
	    s_max_u=max_u_cuts_outb[max_u_cut_lvl]//max_u_med_out
	    s_max_v=max_v_cuts_outb[max_v_cut_lvl]//max_v_med_out
	    s_max_w=max_w_cuts_outb[max_w_cut_lvl]//max_w_med_out
	    
	    heightR1=heightR1_outb
	    heightR2=heightR2_outb
	    heightR3=heightR3_outb
	    
	    radiusR1=radiusR1_outb
	    radiusR2=radiusR2_outb
	    radiusR3=radiusR3_outb

	    min_vz=vz_min_sect_outb
	    max_vz=vz_max_sect_outb

	    p0mean  = p0mean_outb
	    p1mean  = p1mean_outb
	    p2mean  = p2mean_outb
	    p3mean  = p3mean_outb
	    p0sigma = p0sigma_outb
	    p1sigma = p1sigma_outb
	    p_min=0.2381 + 0.11905*ebeam
	    anti_pion_threshold = 0.2

        reg1_min_sec1 = reg1_min_sec1_outb
        reg1_max_sec1 = reg1_max_sec1_outb
        reg2_min_sec1 = reg2_min_sec1_outb
        reg2_max_sec1 = reg2_max_sec1_outb
        reg3_min_sec1 = reg3_min_sec1_outb
        reg3_max_sec1 = reg3_max_sec1_outb
        reg1_min_sec2 = reg1_min_sec2_outb
        reg1_max_sec2 = reg1_max_sec2_outb
        reg2_min_sec2 = reg2_min_sec2_outb
        reg2_max_sec2 = reg2_max_sec2_outb
        reg3_min_sec2 = reg3_min_sec2_outb
        reg3_max_sec2 = reg3_max_sec2_outb
        reg1_min_sec3 = reg1_min_sec3_outb
        reg1_max_sec3 = reg1_max_sec3_outb
        reg2_min_sec3 = reg2_min_sec3_outb
        reg2_max_sec3 = reg2_max_sec3_outb
        reg3_min_sec3 = reg3_min_sec3_outb
        reg3_max_sec3 = reg3_max_sec3_outb
        reg1_min_sec4 = reg1_min_sec4_outb
        reg1_max_sec4 = reg1_max_sec4_outb
        reg2_min_sec4 = reg2_min_sec4_outb
        reg2_max_sec4 = reg2_max_sec4_outb
        reg3_min_sec4 = reg3_min_sec4_outb
        reg3_max_sec4 = reg3_max_sec4_outb
        reg1_min_sec5 = reg1_min_sec5_outb
        reg1_max_sec5 = reg1_max_sec5_outb
        reg2_min_sec5 = reg2_min_sec5_outb
        reg2_max_sec5 = reg2_max_sec5_outb
        reg3_min_sec5 = reg3_min_sec5_outb
        reg3_max_sec5 = reg3_max_sec5_outb
        reg1_min_sec6 = reg1_min_sec6_outb
        reg1_max_sec6 = reg1_max_sec6_outb
        reg2_min_sec6 = reg2_min_sec6_outb
        reg2_max_sec6 = reg2_max_sec6_outb
        reg3_min_sec6 = reg3_min_sec6_outb
        reg3_max_sec6 = reg3_max_sec6_outb

	}
	else if( magnetic_field_config == "inbending" ){
	    println('[ElectronFromEvent::setElectronCutParameters] -> setting parameters for inbending')
	    s_min_u=min_u_cuts_inb[min_u_cut_lvl] //min_u_med_out
	    s_min_v=min_v_cuts_inb[min_v_cut_lvl]
	    s_min_w=min_w_cuts_inb[min_w_cut_lvl]
	    
	    s_max_u=max_u_cuts_inb[max_u_cut_lvl]//max_u_med_out
	    s_max_v=max_v_cuts_inb[max_v_cut_lvl]//max_v_med_out
	    s_max_w=max_w_cuts_inb[max_w_cut_lvl]//max_w_med_out

	    heightR1=heightR1_inb
	    heightR2=heightR2_inb
	    heightR3=heightR3_inb

	    radiusR1=radiusR1_inb
	    radiusR2=radiusR2_inb
	    radiusR3=radiusR3_inb

	    min_vz=vz_min_sect_inb
	    max_vz=vz_max_sect_inb

	    p0mean  = p0mean_inb
	    p1mean  = p1mean_inb
	    p2mean  = p2mean_inb
	    p3mean  = p3mean_inb
	    p0sigma = p0sigma_inb
	    p1sigma = p1sigma_inb	
	    p_min=0.2381 + 0.11905*ebeam
	    anti_pion_threshold=0.2

        reg1_min_sec1 = reg1_min_sec1_inb
        reg1_max_sec1 = reg1_max_sec1_inb
        reg2_min_sec1 = reg2_min_sec1_inb
        reg2_max_sec1 = reg2_max_sec1_inb
        reg3_min_sec1 = reg3_min_sec1_inb
        reg3_max_sec1 = reg3_max_sec1_inb
        reg1_min_sec2 = reg1_min_sec2_inb
        reg1_max_sec2 = reg1_max_sec2_inb
        reg2_min_sec2 = reg2_min_sec2_inb
        reg2_max_sec2 = reg2_max_sec2_inb
        reg3_min_sec2 = reg3_min_sec2_inb
        reg3_max_sec2 = reg3_max_sec2_inb
        reg1_min_sec3 = reg1_min_sec3_inb
        reg1_max_sec3 = reg1_max_sec3_inb
        reg2_min_sec3 = reg2_min_sec3_inb
        reg2_max_sec3 = reg2_max_sec3_inb
        reg3_min_sec3 = reg3_min_sec3_inb
        reg3_max_sec3 = reg3_max_sec3_inb
        reg1_min_sec4 = reg1_min_sec4_inb
        reg1_max_sec4 = reg1_max_sec4_inb
        reg2_min_sec4 = reg2_min_sec4_inb
        reg2_max_sec4 = reg2_max_sec4_inb
        reg3_min_sec4 = reg3_min_sec4_inb
        reg3_max_sec4 = reg3_max_sec4_inb
        reg1_min_sec5 = reg1_min_sec5_inb
        reg1_max_sec5 = reg1_max_sec5_inb
        reg2_min_sec5 = reg2_min_sec5_inb
        reg2_max_sec5 = reg2_max_sec5_inb
        reg3_min_sec5 = reg3_min_sec5_inb
        reg3_max_sec5 = reg3_max_sec5_inb
        reg1_min_sec6 = reg1_min_sec6_inb
        reg1_max_sec6 = reg1_max_sec6_inb
        reg2_min_sec6 = reg2_min_sec6_inb
        reg2_max_sec6 = reg2_max_sec6_inb
        reg3_min_sec6 = reg3_min_sec6_inb
        reg3_max_sec6 = reg3_max_sec6_inb

	}


	println('[ElectronFromEvent::setElectronCutParameters] -> electron pdg pid ' + ebPID )
	println('[ElectronFromEvent::setElectronCutParameters] -> min pcal edep ' + min_ecal_inner_dep )
	println('[ElectronFromEvent::setElectronCutParameters] -> min nphe htcc ' + min_nphe )

	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL U limits per sector : ' + s_min_u )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL U limits per sector : ' + s_max_u )
	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL V limits per sector : ' + s_min_v )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL V limits per sector : ' + s_max_v )
	println('[ElectronFromEvent::setElectronCutParameters] -> min PCAL W limits per sector : ' + s_min_w )
	println('[ElectronFromEvent::setElectronCutParameters] -> max PCAL W limits per sector : ' + s_max_w )

	println('[ElectronFromEvent::setElectronCutParameters] -> DC R1 height ' + heightR1 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R2 height ' + heightR2 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R3 height ' + heightR3 )

	println('[ElectronFromEvent::setElectronCutParameters] -> DC R1 radius ' + radiusR1 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R2 radius ' + radiusR2 )
	println('[ElectronFromEvent::setElectronCutParameters] -> DC R3 radius ' + radiusR3 )

	println('[ElectronFromEvent::setElectronCutParameters] -> min vertexZ limits per sector ' + min_vz )
	println('[ElectronFromEvent::setElectronCutParameters] -> max vertexZ limits per sector ' + max_vz )

	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p0 mean ' + p0mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p1 mean ' + p1mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p2 mean ' + p2mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p3 mean ' + p3mean )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p1 sigma ' + p0sigma )
	println('[ElectronFromEvent::setElectronCutParameters] -> ec sampling fraction p2 sigma ' + p1sigma )

    }

    //////////////////////////////////////////////

    //def heightR2 = 35
    //def radiusR2 = 40

    //def heightR3 = 48
    //def radiusR3 = 49

    ///////////////////////////////////////////////////////////////////////////////////////
        
    def passElectronStatus = { event, index ->
        return (event.status[index] < 0)
    }

    def passElectronEBPIDCut = { event, index ->
	return (event.pid[index] == ebPID)
    }

    def passElectronChargeCut = { event, index ->
	return (event.charge[index] < 0)
    }

    def passElectronNpheCut = { event, index ->
	return (event.cherenkov_status.contains(index)) ? event.nphe[index] > min_nphe : false
    }

    def passElectronVertexCut = { event, index ->
	if (event.pcal_status.contains(index)){	  	
	    def sec = event.pcal_sector[index]-1 // ? do we want the pcal sector or a drift chamber sector?
            return (event.vz[index].with{ it < max_vz[sec] && it > min_vz[sec] })
	}
	return false
    }

    def passElectronTrackQualityCut = { event, index ->
	if(event.cherenkov_status.contains(index) && 
	   event.pcal_status.contains(index) && 
	   event.dc1_status.contains(index) ){
	    if( event.pcal_sector[index] > 0 && event.dc_sector[index] > 0 && event.cherenkov_sector[index] > 0){ //forces tracks to be in Forward detector
		return true	    
	    }
	}
	return false
    }

    def passElectronPCALFiducialCut = { event, index ->
	if (event.pcal_status.contains(index)){	  
	    def pcal_sector = event.pcal_sector[index]-1
	    return (event.pcal_u[index].with{it < s_max_u[pcal_sector] && it > s_min_u[pcal_sector] } &&
		    event.pcal_v[index].with{it < s_max_v[pcal_sector] && it > s_min_v[pcal_sector] } &&
		    event.pcal_w[index].with{it < s_max_w[pcal_sector] && it > s_min_w[pcal_sector] })
	}
	return false
    }

    def passElectronSamplingFractionCut = { event, index ->
	if( event.ecal_inner_status.contains(index) ||  event.ecal_outer_status.contains(index) || event.pcal_status.contains(index) ){
	    def eidep=0
	    def eodep = 0
	    def pcaldep = 0
	    def sector = -1
	    if( event.ecal_inner_status.contains(index) ){
		eidep = event.ecal_inner_energy[index]
		//sector = event.ecal_inner_sector[index] -1
	    }
	    if( event.ecal_outer_status.contains(index) ){
		eodep = event.ecal_outer_energy[index]
		//sector = event.ecal_outer_sector[index]-1
	    }
	    if( event.pcal_status.contains(index) ){
		pcaldep = event.pcal_energy[index]
		sector = event.pcal_sector[index]-1
	    }
	    
	    def edep = eidep + eodep + pcaldep 
	    if( sector >= 0 ){
		
		def p = event.p[index]			
		def mean = p0mean[sector] * (1 + p/Math.sqrt(p*p + p1mean[sector])) + p2mean[sector]*p + p3mean[sector]*p*p
		def sigma = p0sigma[sector] + p1sigma[sector]/Math.sqrt(p)

		def upper_cut = mean + sigma_range * sigma
		def lower_cut = mean - sigma_range * sigma
		
		if( edep/p <= upper_cut && edep/p >= lower_cut ) return true	    
	    }
	}
	return false
    }

    def passElectronAntiPionCut = { event, index ->
	if( event.ecal_inner_status.contains(index) && event.pcal_status.contains(index) ){
	    return  -event.pcal_energy[index]/event.p[index] + anti_pion_threshold < event.ecal_inner_energy[index]/event.p[index]
	}
	return false
    }
	

    def passElectronMinMomentum = { event, index ->
	return (event.p[index] > p_min )
    }

    def passElectronEIEOCut = { event, index ->
	return (event.ecal_inner_status.contains(index)) ? event.ecal_inner_energy[index] > min_ecal_inner_dep : false
    }

    //detector layer r1-12, r2-24, r3-36
    //rotate hit position based on sector
    def rotateDCHitPosition(hit, sec) {
        def ang = Math.toRadians(sec * sect_angle_coverage)
        def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
        def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
        return [x1_rot, y1_rot]
    }

    //define left right 
    def borderDCHitPosition(y_rot, height) {
        def slope = 1 / Math.tan(Math.toRadians(0.5 * sect_angle_coverage))
        def left = (height - slope * y_rot)
        def right = (height + slope * y_rot)
        return [left, right]
    }

    def passElectronDCR1 = { event, index ->
	if (event.dc1_status.contains(index)){
	    def sec = event.dc_sector.get(index)-1
	    def hit = event.dc1.get(index).find{ hit -> hit.layer == 6}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR1[sec])		
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR1 )
            } else {
                return false
            }
	}
	return false
    }

    def passElectronDCR2 = { event, index ->
	if (event.dc2_status.contains(index)){
	    def sec = event.dc_sector[index]-1
	    def hit = event.dc2.get(index).find{ hit -> hit.layer == 18}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR2[sec])
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR2 )
            } else {
                return false
            }
	}
	return false
    }

    def passElectronDCR3 = { event, index ->
	if (event.dc3_status.contains(index)){
	    def sec = event.dc_sector[index]-1
	    def hit = event.dc3.get(index).find{ hit -> hit.layer == 36}
	    if (hit) {
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR3[sec])

                 // from lower line: && x1_rot**2 > radius2_DCr1){ return true }
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) && Math.pow(hit_rotate.get(0),2) > radiusR3 )
            } else {
                return false
            }
	}
	return false
    }

    def passElectronDCR1_updated = { event, index ->
    if (event.dc1_status.contains(index)){
        def sec = event.dc_sector.get(index)-1
        def hit = event.dc1.get(index).find{ hit -> hit.layer == 6}
            if (hit){
                def (theta_DC, phi_DC) = DCthetaphi([hit.x, hit.y, hit.z], sec)
                def (phi_DC_min, phi_DC_max) = DCR1phimaxmin(theta_DC, sec)
                return phi_DC > phi_DC_min && phi_DC < phi_DC_max
            } else {
                return false
            }
    }
    return false
    }

    def passElectronDCR2_updated = { event, index ->
    if (event.dc1_status.contains(index)){
        def sec = event.dc_sector.get(index)-1
        def hit = event.dc2.get(index).find{ hit -> hit.layer == 18}
            if (hit){
                def (theta_DC, phi_DC) = DCthetaphi([hit.x, hit.y, hit.z], sec)
                def (phi_DC_min, phi_DC_max) = DCR2phimaxmin(theta_DC, sec)      
                return phi_DC > phi_DC_min && phi_DC < phi_DC_max
            } else {
                return false
            }
    }
    return false
    }

    def passElectronDCR3_updated = { event, index ->
    if (event.dc1_status.contains(index)){
        def sec = event.dc_sector.get(index)-1
        def hit = event.dc3.get(index).find{ hit -> hit.layer == 36}
            if (hit){
                def (theta_DC, phi_DC) = DCthetaphi([hit.x, hit.y, hit.z], sec)
                def (phi_DC_min, phi_DC_max) = DCR3phimaxmin(theta_DC, sec)      
                return phi_DC > phi_DC_min && phi_DC < phi_DC_max
            } else {
                return false
            }
    }
    return false
    }


    def DCthetaphi(hit, sec){
        def theta_DC = Math.toDegrees(
            Math.acos(hit.get(2)/Math.sqrt(hit.get(0)*hit.get(0) + hit.get(1)*hit.get(1) + hit.get(2)*hit.get(2))))
        def phi_DC_raw = Math.toDegrees(Math.atan2(hit.get(1), hit.get(0)))
        def phi_DC = 0;
        if(sec+1 == 1) phi_DC = phi_DC_raw;
        if(sec+1 == 2) phi_DC = phi_DC_raw - 60;
        if(sec+1 == 3) phi_DC = phi_DC_raw - 120;
        if(sec+1 == 4 && phi_DC_raw > 0) phi_DC = phi_DC_raw - 180;
        if(sec+1 == 4 && phi_DC_raw < 0) phi_DC = phi_DC_raw + 180;
        if(sec+1 == 5) phi_DC = phi_DC_raw + 120;
        if(sec+1 == 6) phi_DC = phi_DC_raw + 60;
        return [theta_DC, phi_DC]
    }

    def DCR1phimaxmin(theta_DC, sec){

    if(sec+1 == 1){
      p0_min = reg1_min_sec1[0]; p1_min = reg1_min_sec1[1]; p2_min = reg1_min_sec1[2]; p3_min = reg1_min_sec1[3]; p4_min = reg1_min_sec1[4];
      p0_max = reg1_max_sec1[0]; p1_max = reg1_max_sec1[1]; p2_max = reg1_max_sec1[2]; p3_max = reg1_max_sec1[3]; p4_max = reg1_max_sec1[4];
    }
    if(sec+1 == 2){
      p0_min = reg1_min_sec2[0]; p1_min = reg1_min_sec2[1]; p2_min = reg1_min_sec2[2]; p3_min = reg1_min_sec2[3]; p4_min = reg1_min_sec2[4];
      p0_max = reg1_max_sec2[0]; p1_max = reg1_max_sec2[1]; p2_max = reg1_max_sec2[2]; p3_max = reg1_max_sec2[3]; p4_max = reg1_max_sec2[4];
    }
    if(sec+1 == 3){
      p0_min = reg1_min_sec3[0]; p1_min = reg1_min_sec3[1]; p2_min = reg1_min_sec3[2]; p3_min = reg1_min_sec3[3]; p4_min = reg1_min_sec3[4];
      p0_max = reg1_max_sec3[0]; p1_max = reg1_max_sec3[1]; p2_max = reg1_max_sec3[2]; p3_max = reg1_max_sec3[3]; p4_max = reg1_max_sec3[4];
    }
    if(sec+1 == 4){
      p0_min = reg1_min_sec4[0]; p1_min = reg1_min_sec4[1]; p2_min = reg1_min_sec4[2]; p3_min = reg1_min_sec4[3]; p4_min = reg1_min_sec4[4];
      p0_max = reg1_max_sec4[0]; p1_max = reg1_max_sec4[1]; p2_max = reg1_max_sec4[2]; p3_max = reg1_max_sec4[3]; p4_max = reg1_max_sec4[4];
    }
    if(sec+1 == 5){
      p0_min = reg1_min_sec5[0]; p1_min = reg1_min_sec5[1]; p2_min = reg1_min_sec5[2]; p3_min = reg1_min_sec5[3]; p4_min = reg1_min_sec5[4];
      p0_max = reg1_max_sec5[0]; p1_max = reg1_max_sec5[1]; p2_max = reg1_max_sec5[2]; p3_max = reg1_max_sec5[3]; p4_max = reg1_max_sec5[4];
    }
    if(sec+1 == 6){
      p0_min = reg1_min_sec6[0]; p1_min = reg1_min_sec6[1]; p2_min = reg1_min_sec6[2]; p3_min = reg1_min_sec6[3]; p4_min = reg1_min_sec6[4];
      p0_max = reg1_max_sec6[0]; p1_max = reg1_max_sec6[1]; p2_max = reg1_max_sec6[2]; p3_max = reg1_max_sec6[3]; p4_max = reg1_max_sec6[4];
    }

    def phi_DC_min = p0_min + p1_min * Math.log(theta_DC) + p2_min * theta_DC + p3_min * theta_DC*theta_DC + p4_min * theta_DC*theta_DC*theta_DC;
    def phi_DC_max = p0_max + p1_max * Math.log(theta_DC) + p2_max * theta_DC + p3_max * theta_DC*theta_DC + p4_max * theta_DC*theta_DC*theta_DC;

    if(phi_DC_min < -25.5) phi_DC_min = -25.5;
    if(phi_DC_max > +25.5) phi_DC_max = +25.5;

    return [phi_DC_min, phi_DC_max];
    }

    def DCR2phimaxmin(theta_DC, sec){

    if(sec+1 == 1){
      p0_min = reg2_min_sec1[0]; p1_min = reg2_min_sec1[1]; p2_min = reg2_min_sec1[2]; p3_min = reg2_min_sec1[3]; p4_min = reg2_min_sec1[4];
      p0_max = reg2_max_sec1[0]; p1_max = reg2_max_sec1[1]; p2_max = reg2_max_sec1[2]; p3_max = reg2_max_sec1[3]; p4_max = reg2_max_sec1[4];
    }
    if(sec+1 == 2){
      p0_min = reg2_min_sec2[0]; p1_min = reg2_min_sec2[1]; p2_min = reg2_min_sec2[2]; p3_min = reg2_min_sec2[3]; p4_min = reg2_min_sec2[4];
      p0_max = reg2_max_sec2[0]; p1_max = reg2_max_sec2[1]; p2_max = reg2_max_sec2[2]; p3_max = reg2_max_sec2[3]; p4_max = reg2_max_sec2[4];
    }
    if(sec+1 == 3){
      p0_min = reg2_min_sec3[0]; p1_min = reg2_min_sec3[1]; p2_min = reg2_min_sec3[2]; p3_min = reg2_min_sec3[3]; p4_min = reg2_min_sec3[4];
      p0_max = reg2_max_sec3[0]; p1_max = reg2_max_sec3[1]; p2_max = reg2_max_sec3[2]; p3_max = reg2_max_sec3[3]; p4_max = reg2_max_sec3[4];
    }
    if(sec+1 == 4){
      p0_min = reg2_min_sec4[0]; p1_min = reg2_min_sec4[1]; p2_min = reg2_min_sec4[2]; p3_min = reg2_min_sec4[3]; p4_min = reg2_min_sec4[4];
      p0_max = reg2_max_sec4[0]; p1_max = reg2_max_sec4[1]; p2_max = reg2_max_sec4[2]; p3_max = reg2_max_sec4[3]; p4_max = reg2_max_sec4[4];
    }
    if(sec+1 == 5){
      p0_min = reg2_min_sec5[0]; p1_min = reg2_min_sec5[1]; p2_min = reg2_min_sec5[2]; p3_min = reg2_min_sec5[3]; p4_min = reg2_min_sec5[4];
      p0_max = reg2_max_sec5[0]; p1_max = reg2_max_sec5[1]; p2_max = reg2_max_sec5[2]; p3_max = reg2_max_sec5[3]; p4_max = reg2_max_sec5[4];
    }
    if(sec+1 == 6){
      p0_min = reg2_min_sec6[0]; p1_min = reg2_min_sec6[1]; p2_min = reg2_min_sec6[2]; p3_min = reg2_min_sec6[3]; p4_min = reg2_min_sec6[4];
      p0_max = reg2_max_sec6[0]; p1_max = reg2_max_sec6[1]; p2_max = reg2_max_sec6[2]; p3_max = reg2_max_sec6[3]; p4_max = reg2_max_sec6[4];
    }

    def phi_DC_min = p0_min + p1_min * Math.log(theta_DC) + p2_min * theta_DC + p3_min * theta_DC*theta_DC + p4_min * theta_DC*theta_DC*theta_DC;
    def phi_DC_max = p0_max + p1_max * Math.log(theta_DC) + p2_max * theta_DC + p3_max * theta_DC*theta_DC + p4_max * theta_DC*theta_DC*theta_DC;

    if(phi_DC_min < -25.5) phi_DC_min = -25.5;
    if(phi_DC_max > +25.5) phi_DC_max = +25.5;

    return [phi_DC_min, phi_DC_max];
    }

    def DCR3phimaxmin(theta_DC, sec){

    if(sec+1 == 1){
      p0_min = reg3_min_sec1[0]; p1_min = reg3_min_sec1[1]; p2_min = reg3_min_sec1[2]; p3_min = reg3_min_sec1[3]; p4_min = reg3_min_sec1[4];
      p0_max = reg3_max_sec1[0]; p1_max = reg3_max_sec1[1]; p2_max = reg3_max_sec1[2]; p3_max = reg3_max_sec1[3]; p4_max = reg3_max_sec1[4];
    }
    if(sec+1 == 2){
      p0_min = reg3_min_sec2[0]; p1_min = reg3_min_sec2[1]; p2_min = reg3_min_sec2[2]; p3_min = reg3_min_sec2[3]; p4_min = reg3_min_sec2[4];
      p0_max = reg3_max_sec2[0]; p1_max = reg3_max_sec2[1]; p2_max = reg3_max_sec2[2]; p3_max = reg3_max_sec2[3]; p4_max = reg3_max_sec2[4];
    }
    if(sec+1 == 3){
      p0_min = reg3_min_sec3[0]; p1_min = reg3_min_sec3[1]; p2_min = reg3_min_sec3[2]; p3_min = reg3_min_sec3[3]; p4_min = reg3_min_sec3[4];
      p0_max = reg3_max_sec3[0]; p1_max = reg3_max_sec3[1]; p2_max = reg3_max_sec3[2]; p3_max = reg3_max_sec3[3]; p4_max = reg3_max_sec3[4];
    }
    if(sec+1 == 4){
      p0_min = reg3_min_sec4[0]; p1_min = reg3_min_sec4[1]; p2_min = reg3_min_sec4[2]; p3_min = reg3_min_sec4[3]; p4_min = reg3_min_sec4[4];
      p0_max = reg3_max_sec4[0]; p1_max = reg3_max_sec4[1]; p2_max = reg3_max_sec4[2]; p3_max = reg3_max_sec4[3]; p4_max = reg3_max_sec4[4];
    }
    if(sec+1 == 5){
      p0_min = reg3_min_sec5[0]; p1_min = reg3_min_sec5[1]; p2_min = reg3_min_sec5[2]; p3_min = reg3_min_sec5[3]; p4_min = reg3_min_sec5[4];
      p0_max = reg3_max_sec5[0]; p1_max = reg3_max_sec5[1]; p2_max = reg3_max_sec5[2]; p3_max = reg3_max_sec5[3]; p4_max = reg3_max_sec5[4];
    }
    if(sec+1 == 6){
      p0_min = reg3_min_sec6[0]; p1_min = reg3_min_sec6[1]; p2_min = reg3_min_sec6[2]; p3_min = reg3_min_sec6[3]; p4_min = reg3_min_sec6[4];
      p0_max = reg3_max_sec6[0]; p1_max = reg3_max_sec6[1]; p2_max = reg3_max_sec6[2]; p3_max = reg3_max_sec6[3]; p4_max = reg3_max_sec6[4];
    }

    def phi_DC_min = p0_min + p1_min * Math.log(theta_DC) + p2_min * theta_DC + p3_min * theta_DC*theta_DC + p4_min * theta_DC*theta_DC*theta_DC;
    def phi_DC_max = p0_max + p1_max * Math.log(theta_DC) + p2_max * theta_DC + p3_max * theta_DC*theta_DC + p4_max * theta_DC*theta_DC*theta_DC;

    if(phi_DC_min < -25.5) phi_DC_min = -25.5;
    if(phi_DC_max > +25.5) phi_DC_max = +25.5;

    return [phi_DC_min, phi_DC_max];
    }


}
