package sangbaek
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class draw_dcr {
	def hists = new ConcurrentHashMap()
	def h_DCR_hits1 = {new H2F("$it","$it",200,-300,300, 200,-300,300)}
	def h_DCR_hits2 = {new H2F("$it","$it",200,-400,400, 200,-400,400)}
	def h_DCR_hits3 = {new H2F("$it","$it",200,-500,500, 200,-500,500)}

	def processEvent(event) {
		(0..<event.npart).each{ index ->
            if (event.dc1_status.contains(index)){
                    def sec = event.dc_sector.get(index)-1
                    def hit1 = event.dc1.get(index)
                    if (hit1){
                            hit1.each{ hists.computeIfAbsent("/DCR1/layer"+it.layer+"/sec_"+(sec+1),h_DCR_hits1).fill(it.x, it.y)}
                    }   
            }   

            if (event.dc2_status.contains(index)){
                    def sec = event.dc_sector.get(index)-1
                    def hit2 = event.dc2.get(index)
                    if (hit2){
                            hit2.each{ hists.computeIfAbsent("/DCR2/layer"+it.layer+"/sec_"+(sec+1),h_DCR_hits2).fill(it.x, it.y)}
                    }   
            }   

            if (event.dc3_status.contains(index)){
                    def sec = event.dc_sector.get(index)-1
                    def hit3 = event.dc3.get(index)
                    if (hit3){
                            hit3.each{ hists.computeIfAbsent("/DCR3/layer"+it.layer+"/sec_"+(sec+1),h_DCR_hits3).fill(it.x, it.y)}
                    }   
            }   

		}
	}
}