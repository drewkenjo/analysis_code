package pid.sangbaek

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import pid.electron.ElectronFromEvent
import pid.electron.ElectronSelector
import event.Event
import org.jlab.clas.physics.Vector3

class electron{

  def event

  def electron_candidate = new ElectronFromEvent()
  def electron_selector = new ElectronSelector()
  def electronCutStrategies_Brandon
  def electronCutStrategies_Custom

  def electronCutResults_Brandon
  def electronCutResults_Custom

  def electron(){
    this.initalizeCustomElecCuts()
  }

  def electron(event){
    this.event = event
    this.initalizeCustomElecCuts()
    this.getGoodElectron(event)
    this.getGoodElectronCustom(event)
  }

  def applyCuts_Brandon(event){
    this.getGoodElectron(event)
    return this.electronCutResults_Brandon
  }

  def applyCuts_Custom(event){
    this.getGoodElectronCustom(event)
    return this.electronCutResults_Custom
  }


  def initalizeCustomElecCuts(){
    this.electronCutStrategies_Brandon = [
      this.electron_candidate.passElectronStatus,
      this.electron_candidate.passElectronChargeCut,
      this.electron_candidate.passElectronTrackQualityCut,
      this.electron_candidate.passElectronMinMomentum,
      this.electron_candidate.passElectronEBPIDCut,
      this.electron_candidate.passElectronSamplingFractionCut,
      this.electron_candidate.passElectronNpheCut,
      this.electron_candidate.passElectronVertexCut,
      this.electron_candidate.passElectronPCALFiducialCut,
      this.electron_candidate.passElectronEIEOCut,
      this.electron_candidate.passElectronDCR1,
      this.electron_candidate.passElectronDCR2,
      this.electron_candidate.passElectronDCR3,
      this.electron_candidate.passElectronAntiPionCut
    ]

    this.electronCutStrategies_Custom = [
      this.find_byFTOF,
      this.find_byMOM
    ]

    def field_setting = "inbending"
    // cut lvl meanings: 0 loose, 1 med, 2 tight
    def el_cut_strictness_lvl=["ecal_cut_lvl":1,
               "nphe_cut_lvl":1,
               "vz_cut_lvl":1,
               "min_u_cut_lvl":1,
               "min_v_cut_lvl":1,
               "min_w_cut_lvl":1,
               "max_u_cut_lvl":1,
               "max_v_cut_lvl":1,
               "max_w_cut_lvl":1,
               "dcr1_cut_lvl":1,
               "dcr2_cut_lvl":1,
               "dcr3_cut_lvl":1,
               "anti_pion_cut_lvl":1
    ]
    this.electron_selector.setElectronCutStrictness(el_cut_strictness_lvl)
    this.electron_selector.setCutParameterFromMagField("inbending")

    this.electron_candidate.setElectronCutStrictness(el_cut_strictness_lvl)
    this.electron_candidate.setElectronCutParameters("inbending")

  }


  def getGoodElectron(event){
    //return a list of REC::Particle indices for tracks passing all electron cuts
    def el_cut_result = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, this.electronCutStrategies_Brandon.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
    this.electronCutResults_Brandon = el_cut_result.findResults{el_indx, cut_result -> !cut_result.contains(false) ? el_indx : null}
  }

  def getGoodElectronCustom(event){
    this.electronCutResults_Custom = this.electronCutResults_Brandon.findResults{ index ->
      this.electronCutStrategies_Custom.collect{custom_test -> !custom_test(event,index)}.contains(false)? index : null
    }
  }

  // FTOF Hit Response
  def find_byFTOF = {event, index ->
    return event.tof_status.contains(index)
  }

  // momentum cut
  def find_byMOM = {event, index ->
    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
    def p = lv.mag()
    def vz = event.vz[index]
    def theta = Math.toDegrees(lv.theta())
    def phi = Math.toDegrees(lv.phi())
    return p > 1.5 && theta>17*(1-p/7)
  }
}