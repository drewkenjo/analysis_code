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
import pid.sangbaek.electron

def field_setting = "inbending"
// cut lvl meanings: 0 loose, 1 med, 2 tight
el_cut_strictness_lvl=["ecal_cut_lvl":1,
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

def hist_ele_theta_mom = [:].withDefault{new H2F("hist_ele_theta_mom_$it", "momentum vs theta", 250,0,30, 250,0,11)}

//use these two lines for the second and third method
def ele_selector = new ElectronSelector()
//if you want to use the ElectronSelector class
//ele_selector.initializeCuts()
//use the two methods below
//must be called in this order
ele_selector.setElectronCutStrictness(el_cut_strictness_lvl)
ele_selector.setCutParameterFromMagField("inbending")

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

def electron_ind = new electron()
   
for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    while(reader.hasEvent()) {
		def data_event = reader.getNextEvent()
		def event = EventConverter.convert(data_event)

		def good_el_with_cuts_Brandon = electron_ind.applyCuts_Brandon(event)
		def good_el_with_cuts_Custom = electron_ind.applyCuts_Custom(event)

		good_el_with_cuts_Brandon.each{eleind->
	        def ele = new Particle(11, *[event.px, event.py, event.pz].collect{it[eleind]})
			def ele_p = ele.p()
			def ele_theta = Math.toDegrees(ele.theta())
			hist_ele_theta_mom['Brandon'+event.pcal_sector[eleind]].fill(ele_theta,ele_p)
		}

		good_el_with_cuts_Custom.each{ eleind ->
	        def ele = new Particle(11, *[event.px, event.py, event.pz].collect{it[eleind]})
			def ele_p = ele.p()
			def ele_theta = Math.toDegrees(ele.theta())
			hist_ele_theta_mom['Custom'+event.pcal_sector[eleind]].fill(ele_theta,ele_p)
		}
    }
    reader.close()
}

def out = new TDirectory()
out.mkdir('/electron')
out.cd('/electron')
hist_ele_theta_mom.values().each{out.addDataSet(it)}
out.writeFile('electron_out.hipo')