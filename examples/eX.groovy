import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import utils.KinTool
import pid.electron.ElectronSelector
import pid.proton.ProtonFromEvent
import pid.proton.ProtonSelector
import pid.gamma.GammaFromEvent
import pid.gamma.GammaSelector

def field_setting = "outbending"
// cut lvl meanings: 0 loose, 1 med, 2 tight
el_cut_strictness_lvl=["ecal_cut_lvl":1,
		       "nphe_cut_lvl":1,
		       "vz_cut_lvl":1,
		       "anti_pion_cut_lvl":1
]


//use these two lines for the second and third method
def ele_selector = new ElectronSelector()
//if you want to use the ElectronSelector class
//ele_selector.initializeCuts()
//use the two methods below
//must be called in this order
ele_selector.setElectronCutStrictness(el_cut_strictness_lvl)
ele_selector.setCutParameterFromMagField(field_setting)

//use two lines below for first method
def electron = new ElectronFromEvent();
//if you want to do it manually use the two lines below
electron.setElectronCutStrictness(el_cut_strictness_lvl)
electron.setElectronCutParameters(field_setting)
 
def myElectronCutStrategies = [
    electron.passElectronStatus,
    electron.passElectronChargeCut,
    electron.passElectronTrackQualityCut,
    electron.passElectronMinMomentum,
    electron.passElectronEBPIDCut,
    electron.passElectronSamplingFractionCut,
    electron.passElectronNpheCut,
    electron.passElectronVertexCut,
    electron.passElectronPCALFiducialCut,
    electron.passElectronEIEOCut,
    electron.passElectronDCR1,
    electron.passElectronDCR2,
    electron.passElectronDCR3,
    electron.passElectronAntiPionCut
]

//Similary, we can keep the same structures for proton
//use two lines below for first method
//there is no cut level for proton at this stage
def pro_selector = new ProtonSelector()
def proton = new ProtonFromEvent();
def myProtonCutStrategies = [
	proton.passProtonEBPIDCut,
	proton.passProtonDCR1,
	proton.passProtonDCR2,
	proton.passProtonDCR3,
    proton.passProtonTrackQuality,
    proton.passProtonCDPolarAngleCut,
    proton.passProtonVertexCut
]
pro_selector.setCutParameterFromMagField(field_setting)
proton.setProtonCutParameters(field_setting)

//the same for gamma
//use these two lines for the second and third method
def gam_selector = new GammaSelector()
//use two lines below for first method
def gamma = new GammaFromEvent();
 
def myGammaCutStrategies = [
    gamma.passGammaEBPIDCut,
    gamma.passGammaPCALFiducialCut,
    gamma.passGammaBetaCut
]
   
for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    while(reader.hasEvent()) {
	def data_event = reader.getNextEvent()
	def event = EventConverter.convert(data_event)

	// first method for selecting electrons illustrates what the ElectronSelector class is doing under the hood.
	def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  
	my_el_cuts.each{ index, value ->
		if (!value.contains(false)){
		    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		    def p = lv.mag()
		    def vz = event.vz[index]
		    def theta = Math.toDegrees(lv.theta())
		    def phi = Math.toDegrees(lv.phi())
	    	//do other stuff here
	    }
	}

	// second method here will return a list of indicies for tracks passing all electron cuts.
	// decide which of the electron candidates you want here
	def my_good_el = ele_selector.getGoodElectron(event)	
	// third method will return a map with key as the REC::Particle index and value as the list of booleans describing if the track passed the cut or not.
	def my_good_el_with_cuts = ele_selector.getGoodElectronWithCuts(event)

	// // do the same for proton
	def my_pro_cuts = (0..<event.npart).findAll{event.charge[it]>0}.collect{ ii -> [ii, myProtonCutStrategies.collect{ pro_test -> pro_test(event,ii) } ] }.collectEntries()	  
	my_pro_cuts.each{ index, value ->
		if (!value.contains(false)){
		    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		    def p = lv.mag()
		    def vz = event.vz[index]
		    def theta = Math.toDegrees(lv.theta())
		    def phi = Math.toDegrees(lv.phi())
		    //do other stuff here
		}
	}

	def my_good_pro = pro_selector.getGoodProton(event)	
	def my_good_pro_with_cuts = pro_selector.getGoodProtonWithCuts(event)

	// do the same for gamma
	def my_gam_cuts = (0..<event.npart).findAll{event.charge[it]==0}.collect{ ii -> [ii, myGammaCutStrategies.collect{ gam_test -> gam_test(event,ii) } ] }.collectEntries()	  
	my_gam_cuts.each{ index, value ->
		if (!value.contains(false)){
		    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		    def p = lv.mag()
		    def vz = event.vz[index]
		    def theta = Math.toDegrees(lv.theta())
		    def phi = Math.toDegrees(lv.phi())
		    //do other stuff here
		}
	}

	def my_good_gam = gam_selector.getGoodGamma(event)	
	def my_good_gam_with_cuts = gam_selector.getGoodGammaWithCuts(event)
    }


}

reader.close()
