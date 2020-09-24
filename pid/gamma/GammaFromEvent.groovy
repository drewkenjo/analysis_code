package pid.gamma

import org.jlab.detector.base.DetectorType
import event.Event

class GammaFromEvent {

  def ebeam = 10.6
  def ebPID = 22

  def min_v = 9
  def min_w = 9

  def passGammaEBPIDCut = { event, index ->
    return (event.pid[index] == ebPID)
  }

  def passGammaPCALFiducialCut = { event, index ->

    if (event.status[index]<2000) return true

    if (event.pcal_status.contains(index)){   
        return ( event.pcal_v[index] > min_v  && event.pcal_w[index] > min_w)
    }
    return false
  }

  def passGammaBetaCut = { event, index ->
    event.beta[index] < 1.1 && event.beta[index]>0.9
  }

}