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

def hist_ele_theta_mom = [:].withDefault{new H2F("hist_ele_theta_mom_$it", "momentum vs theta", 250,0,30, 250,0,11)}
 
def electron_ind = new electron()
   
for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    while(reader.hasEvent()) {
		def data_event = reader.getNextEvent()
		def event = EventConverter.convert(data_event)

		def good_el_with_cuts_Brandon = electron_ind.applyCuts_Brandon(event)
		def good_el_with_cuts_Custom = electron_ind.applyCuts_Custom(event)

		if(good_el_with_cuts_Brandon){
			println(event.event_number)
			println(good_el_with_cuts_Brandon)
		}
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