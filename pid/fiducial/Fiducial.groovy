// Translation of fiducial cuts from Stefan Diehl (sdielhl@jlab.org)
// , cxx to groovy,
// by Sangbaek Lee (sangbaek@jlab.org)
// 
package pid.fiducial
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

class Fiducial {

// define global 


/// PCAL fiducial cuts:

def EC_hit_position_fiducial_cut = {event, index, polarity, mode ->


  // choose polarity
  def inbending = polarity
  def outbending = !polarity

  ///////////////////////////
  // choose tightness of cut, arg = mode
  def tight = false
  def medium = false
  def loose = false

  if (mode == 'tight')  tight = true
  else if (mode == 'medium') medium = true
  else if (mode == 'loose') loose = true
  else{
  	println('Choose at least one fiducial cut mode correctly among tight, medium, and loose. Halt!')
  	return
  }
  //////////////////////////

// Cut using the natural directions of the scintillator bars/ fibers:

  def u = event.pcal_u[index]
  def v = event.pcal_v[index]
  def w = event.pcal_w[index]
   
  /// v + w is going from the side to the back end of the PCAL, u is going from side to side
  /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.

  ///////////////////////////////////////////////////////////////////
  /// inbending:

  // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
  def min_u_tight_inb = [42, 32, 38, 27.5, 32, 29]
  def min_u_med_inb   = [33, 26, 34, 22,   27, 25]
  def min_u_loose_inb = [28, 22, 30, 18,   22, 22]

  // u shows a slight fall of the sampling fraction for high values
  def max_u_tight_inb = [398, 398, 398, 398, 398, 398] 
  def max_u_med_inb   = [408, 408, 408, 408, 408, 408] 
  def max_u_loose_inb = [420, 420, 420, 420, 420, 420] 

  // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
  def min_v_tight_inb = [18.0, 12.0, 19.5,  15.5,  20.0, 13.0]
  def min_v_med_inb   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
  def min_v_loose_inb = [10.25, 8.0, 12.75, 12.5,  13.25, 9.0]

  // the maximum of v is never reached
  def max_v_tight_inb = [400, 400, 400, 400, 400, 400]
  def max_v_med_inb   = [400, 400, 400, 400, 400, 400]
  def max_v_loose_inb = [400, 400, 400, 400, 400, 400]

  // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
  def min_w_tight_inb = [14.0, 18.7, 18.7,  12.0, 16.0, 13.0]
  def min_w_med_inb   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
  def min_w_loose_inb = [7.25, 11.0, 13.0,  5.5,  10.0, 6.0]

  // the maximum of w is never reached
  def max_w_tight_inb = [400, 400, 400, 400, 400, 400]
  def max_w_med_inb   = [400, 400, 400, 400, 400, 400]
  def max_w_loose_inb = [400, 400, 400, 400, 400, 400]


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

  // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
  def min_v_tight_out = [18.0, 12.0, 19.5,  15.5,  20.0, 13.0]
  def min_v_med_out   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
  def min_v_loose_out = [10.25, 8.0, 12.75, 12.5,  13.25, 9.0]

  // the maximum of v is never reached
  def max_v_tight_out = [400, 400, 400, 400, 400, 400]
  def max_v_med_out   = [400, 400, 400, 400, 400, 400]
  def max_v_loose_out = [400, 400, 400, 400, 400, 400]

  // tight: only background outside 3 sigma, medium: 10 % outside 3 sigma, loose: 50% outside 3 sigma 
  def min_w_tight_out = [14.0, 18.7, 18.7,  12.0, 16.0, 13.0]
  def min_w_med_out   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
  def min_w_loose_out = [7.25, 11.0, 13.0,  5.5,  10.0, 6.0]

  // the maximum of w is never reached
  def max_w_tight_out = [400, 400, 400, 400, 400, 400]
  def max_w_med_out   = [400, 400, 400, 400, 400, 400]
  def max_w_loose_out = [400, 400, 400, 400, 400, 400]

  //////////////////////////////////////////////////////////////

  def min_u = 0
  def max_u = 0
  def min_v = 0
  def max_v = 0
  def min_w = 0
  def max_w = 0


  (0..<6).each{k->  
    if(event.pcal_sector[index]-1 == k && inbending == true){
      if(tight == true){
        min_u = min_u_tight_inb[k]
        max_u = max_u_tight_inb[k]
        min_v = min_v_tight_inb[k]
        max_v = max_v_tight_inb[k]
        min_w = min_w_tight_inb[k]
        max_w = max_w_tight_inb[k]
      }
      if(medium == true){
        min_u = min_u_med_inb[k]
        max_u = max_u_med_inb[k]
        min_v = min_v_med_inb[k]
        max_v = max_v_med_inb[k]
        min_w = min_w_med_inb[k]
        max_w = max_w_med_inb[k]
      }
      if(loose == true){
        min_u = min_u_loose_inb[k]; max_u = max_u_loose_inb[k];
        min_v = min_v_loose_inb[k]; max_v = max_v_loose_inb[k];
        min_w = min_w_loose_inb[k]; max_w = max_w_loose_inb[k];
      }
    }
    if(event.pcal_sector[index]-1 == k && outbending == true){
      if(tight == true){
        min_u = min_u_tight_out[k]; max_u = max_u_tight_out[k];
        min_v = min_v_tight_out[k]; max_v = max_v_tight_out[k];
        min_w = min_w_tight_out[k]; max_w = max_w_tight_out[k];
      }
      if(medium == true){
        min_u = min_u_med_out[k]; max_u = max_u_med_out[k];
        min_v = min_v_med_out[k]; max_v = max_v_med_out[k];
        min_w = min_w_med_out[k]; max_w = max_w_med_out[k];
      }
      if(loose == true){
        min_u = min_u_loose_out[k]; max_u = max_u_loose_out[k];
        min_v = min_v_loose_out[k]; max_v = max_v_loose_out[k];
        min_w = min_w_loose_out[k]; max_w = max_w_loose_out[k];
      }
    }
  }

  if(v > min_v && v < max_v && w > min_w && w < max_w) return true;
  else return false;
}