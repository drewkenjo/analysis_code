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

  def initalizeCustomElecCuts(){
    this.electronCutStrategies_Brandon = [
      electronCuts.passElectronStatus,
      electronCuts.passElectronChargeCut,
      electronCuts.passElectronTrackQualityCut,
      electronCuts.passElectronMinMomentum,
      electronCuts.passElectronEBPIDCut,
      electronCuts.passElectronSamplingFractionCut,
      electronCuts.passElectronNpheCut,
      electronCuts.passElectronVertexCut,
      electronCuts.passElectronPCALFiducialCut,
      electronCuts.passElectronEIEOCut,
      electronCuts.passElectronDCR1,
      electronCuts.passElectronDCR2,
      electronCuts.passElectronDCR3,
      electronCuts.passElectronAntiPionCut
    ]

    this.electronCutStrategies_Custom = [
      this.find_byFTOF,
      this.find_byMOM
    ]
  }


  def getGoodElectron(event){
    //return a list of REC::Particle indices for tracks passing all electron cuts
    def el_cut_result = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, electronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
    this.electronCutResults_Brandon = el_cut_result.findResults{el_indx, cut_result -> !cut_result.contains(false) ? el_indx : null}
  }

  def getGoodElectronCustom(event){
    this.electronCutResults_Custom = this.electronCutResults_Brandon.findResults{ index ->
      this.electronCutStrategies_Custom.collect{custom_test -> custom_test(event,index)}.contains(false)? index : null
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
    return mom > 1.5 && theta>17*(1-mom/7)
  }
}