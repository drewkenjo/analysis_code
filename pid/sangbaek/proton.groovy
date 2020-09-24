package pid.sangbaek

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import pid.proton.ProtonFromEvent
import event.Event
import org.jlab.clas.physics.Vector3

class proton{

  def event
  def proton_candidate = new ProtonFromEvent()

  def protonCutStrategies_Stefan
  def protonCutStrategies_Custom

  def protonCutResults_Stefan
  def protonCutResults_Custom

  def proton(){
    this.initalizeCustomProCuts()
  }

  def proton(event){
    this.event = event
    this.initalizeCustomProCuts()
    this.getGoodProton(event)
    this.getGoodProtonCustom(event)
  }

  def applyCuts_Stefan(event){
    this.getGoodProton(event)
    return this.protonCutResults_Stefan
  }

  def applyCuts_Custom(event){
    this.getGoodProtonCustom(event)
    return this.protonCutResults_Custom
  }


  def initalizeCustomProCuts(){
    this.protonCutStrategies_Stefan = [
      this.proton_candidate.passProtonEBPIDCut,
      this.proton_candidate.passProtonDCR1,
      this.proton_candidate.passProtonDCR2,
      this.proton_candidate.passProtonDCR3,
      this.proton_candidate.passProtonTrackQuality,
      this.proton_candidate.passProtonCDPolarAngleCut,
      this.proton_candidate.passProtonVertexCut
    ]

    // this.protonCutStrategies_Custom = [
    //   this.find_byFTOF,
    //   this.find_byMOM
    // ]

    // def field_setting = "inbending"
    // cut lvl meanings: 0 loose, 1 med, 2 tight
    // def pro_cut_strictness_lvl=[
    //            "dcr1_cut_lvl":1,
    //            "dcr2_cut_lvl":1,
    //            "dcr3_cut_lvl":1,
    // ]

    // this.proton_candidate.setProtonCutStrictness(pro_cut_strictness_lvl)
    // this.proton_candidate.setProtonCutParameters("inbending")

  }


  def getGoodProton(event){
    //return a list of REC::Particle indices for tracks passing all proton cuts
    def pro_cut_result = (0..<event.npart).findAll{event.charge[it]>0}.collect{ ii -> [ii, this.protonCutStrategies_Stefan.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
    this.protonCutResults_Stefan = pro_cut_result.findResults{el_indx, cut_result -> !cut_result.contains(false) ? el_indx : null}
  }

  def getGoodProtonCustom(event){
    this.protonCutResults_Custom = this.protonCutResults_Stefan.findResults{ index ->
      this.protonCutStrategies_Custom.collect{custom_test -> !custom_test(event,index)}.contains(false)? index : null
    }
  }

  // // FTOF Hit Response
  // def find_byFTOF = {event, index ->
  //   return event.tof_status.contains(index)
  // }

  // // momentum cut
  // def find_byMOM = {event, index ->
  //   def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
  //   def p = lv.mag()
  //   def vz = event.vz[index]
  //   def theta = Math.toDegrees(lv.theta())
  //   def phi = Math.toDegrees(lv.phi())
  //   return p > 1.5 && theta>17*(1-p/7)
  // }
}