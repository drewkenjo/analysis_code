package pid.proton

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import pid.proton.ProtonFromEvent
import event.Event
import org.jlab.clas.physics.Vector3

class ProtonSelector{

  def event
  def protonCuts = new ProtonFromEvent()
  def protonCutStrategies

  def ProtonSelector(){
    this.initalizeCustomProCuts()
  }

  def initalizeCustomProCuts(){
    this.protonCutStrategies = [
      this.protonCuts.passProtonEBPIDCut,
      this.protonCuts.passProtonDCR1,
      this.protonCuts.passProtonDCR2,
      this.protonCuts.passProtonDCR3
    ]
  }

  def getGoodProtonWithCuts(event){
    //return a list of REC::Particle indices for tracks passing all proton cuts
    (0..<event.npart).findAll{event.charge[it]>0}.collect{ ii -> [ii, this.protonCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
  }

  def getGoodProton(event){
    //return a list of REC::Particle indices for tracks passing all proton cuts
    def pro_cut_result = (0..<event.npart).findAll{event.charge[it]>0}.collect{ ii -> [ii, this.protonCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
    return pro_cut_result.findResults{el_indx, cut_result -> !cut_result.contains(false) ? el_indx : null}
  }
}