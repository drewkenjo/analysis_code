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

    def min_u_tight_inb = [42, 32, 38, 27.5, 32, 29]
    def min_u_med_inb   = [33, 26, 34, 22,   27, 25]
    def min_u_loose_inb = [28, 22, 30, 18,   22, 22]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_inb = [398, 398, 398, 398, 398, 398]
    def max_u_med_inb   = [408, 408, 408, 408, 408, 408]
    def max_u_loose_inb = [420, 420, 420, 420, 420, 420] 
    
    def min_u_cuts_inb = [min_u_loose_inb, min_u_med_inb, min_u_tight_inb]
    def max_u_cuts_inb = [max_u_loose_inb, max_u_med_inb, max_u_tight_inb]

    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_inb = [18.0, 12.0, 19.5,  15.5,  20.0, 13.0]
    def min_v_med_inb   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
    def min_v_loose_inb = [10.25, 8.0, 12.75, 12.5,  13.25, 9.0]
    
    // the maximum of v is never reached
    def max_v_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_v_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_inb = [min_v_loose_inb, min_v_med_inb, min_v_tight_inb]
    def max_v_cuts_inb = [max_v_loose_inb, max_v_med_inb, max_v_tight_inb]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_w_tight_inb = [14.0, 18.7, 18.7,  12.0, 16.0, 13.0]
    def min_w_med_inb   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    def min_w_loose_inb = [7.25, 11.0, 13.0,  5.5,  10.0, 6.0]
    
    // the maximum of w is never reached
    def max_w_tight_inb = [400, 400, 400, 400, 400, 400]
    def max_w_med_inb   = [400, 400, 400, 400, 400, 400]
    def max_w_loose_inb = [400, 400, 400, 400, 400, 400]

    def min_w_cuts_inb = [min_w_loose_inb, min_w_med_inb, min_w_tight_inb]
    def max_w_cuts_inb = [max_w_loose_inb, max_w_med_inb, max_w_tight_inb]
   
    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_u_tight_out = [42, 32, 38, 27.5, 32, 29]
    def min_u_med_out   = [33, 26, 34, 22,   27, 25]
    def min_u_loose_out = [28, 22, 30, 18,   22, 22]
    
    // u shows a slight fall of the sampling fraction for high values
    def max_u_tight_out = [398, 398, 398, 398, 398, 398]
    def max_u_med_out   = [408, 408, 408, 408, 408, 408] 
    def max_u_loose_out = [420, 420, 420, 420, 420, 420] 

    def min_u_cuts_outb = [min_u_loose_out, min_u_med_out, min_u_tight_out]
    def max_u_cuts_outb = [max_u_loose_out, max_u_med_out, max_u_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_v_tight_out = [18.0, 12.0, 19.5,  15.5,  20.0, 13.0]
    def min_v_med_out   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
    def min_v_loose_out = [10.25, 8.0, 12.75, 12.5,  13.25, 9.0]
    
    // the maximum of v is never reached
    def max_v_tight_out = [400, 400, 400, 400, 400, 400]
    def max_v_med_out   = [400, 400, 400, 400, 400, 400]
    def max_v_loose_out = [400, 400, 400, 400, 400, 400]

    def min_v_cuts_outb = [min_v_loose_out, min_v_med_out, min_v_tight_out]
    def max_v_cuts_outb = [max_v_loose_out, max_v_med_out, max_v_tight_out]
    
    // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
    def min_w_tight_out = [14.0, 18.7, 18.7,  12.0, 16.0, 13.0]
    def min_w_med_out   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    def min_w_loose_out = [7.25, 11.0, 13.0,  5.5,  10.0, 6.0]
    
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
	}
	else if( dcr1_cut_lvl == 1 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_med }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_med }
	    radiusR1_inb+=dcR1_med
	    radiusR1_outb+=dcR1_med
	}
	else if( dcr1_cut_lvl == 2 ) {
	    heightR1_inb = heightR1_inb.collect{it + dcR1_tight }
	    heightR1_outb = heightR1_outb.collect{it + dcR1_tight }
	    radiusR1_inb+=dcR1_tight
	    radiusR1_outb+=dcR1_tight
	}
	if( dcr2_cut_lvl == 0 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_loose }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_loose }
	    radiusR2_inb+=dcR2_loose
	    radiusR2_outb+=dcR2_loose
	}
	else if( dcr2_cut_lvl == 1 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_med }
	    heightR2_outb = heightR2_outb.collect{it + dcR2_med }
	    radiusR2_inb+=dcR2_med
	    radiusR2_outb+=dcR2_med
	}
	else if( dcr2_cut_lvl == 2 ) {
	    heightR2_inb = heightR2_inb.collect{it + dcR2_tight }
	     heightR2_outb = heightR2_outb.collect{it + dcR2_tight }
	    radiusR2_inb+=dcR2_tight
	    radiusR2_outb+=dcR2_tight
	}

	if( dcr3_cut_lvl == 0 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_loose }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_loose }
	    radiusR3_inb+=dcR3_loose
	    radiusR3_outb+=dcR3_loose
	}
	else if( dcr3_cut_lvl == 1 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_med }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_med }
	    radiusR3_inb+=dcR3_med
	    radiusR3_outb+=dcR3_med
	}
	else if( dcr3_cut_lvl == 2 ) {
	    heightR3_inb = heightR3_inb.collect{it + dcR3_tight }
	    heightR3_outb = heightR3_outb.collect{it + dcR3_tight }
	    radiusR3_inb+=dcR3_tight
	    radiusR3_outb+=dcR3_tight
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
	    if( event.pcal_sector[index]-1 > 0 && event.dc_sector[index]-1 > 0 && event.cherenkov_sector[index]-1 > 0){ //forces tracks to be in Forward detector
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
	    def hit = event.dc1.get(index).find{ hit -> hit.layer == 12}
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
	    def hit = event.dc2.get(index).find{ hit -> hit.layer == 24}
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

    

}
