package sangbaek
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class draw_dcr {
	def hists = new ConcurrentHashMap()
	def h_DCR_hits = {new H2F("$it","$it",200,-600,600, 200,-600,600)}

	def processEvent(event) {
		(0..<event.npart).each{ index ->
			if (event.dc1_status.contains(index)){
				def sec = event.dc_sector.get(index)-1
				def hit1 = event.dc1.get(index).find{ hit -> hit.layer == 12}
				if (hit1) hists.computeIfAbsent("/DCR1/sec_"+(sec+1),h_DCR_hits).fill(hit1.x, hit1.y)
			}

			if (event.dc2_status.contains(index)){
				def sec = event.dc_sector.get(index)-1
				def hit2 = event.dc2.get(index).find{ hit -> hit.layer == 24}
				if (hit2) hists.computeIfAbsent("/DCR2/sec_"+(sec+1),h_DCR_hits).fill(hit2.x, hit2.y)
			}

			if (event.dc1_status.contains(index)){
				def sec = event.dc_sector.get(index)-1
				def hit3 = event.dc3.get(index).find{ hit -> hit.layer == 36}
				if (hit3) hists.computeIfAbsent("/DCR3/sec_"+(sec+1),h_DCR_hits).fill(hit3.x, hit3.y)
			}
		}
	}
}