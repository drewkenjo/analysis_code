package pid.gamma
import event.Event
import pid.gamma.GammaFromEvent

class GammaSelector {

    def gammaCuts = new GammaFromEvent()
    def gammaCutStrategies

    GammaSelector(){
	this.initializeCuts()
    }
    
    def initializeCuts(){
 	this.gammaCutStrategies = [
	    gammaCuts.passGammaEBPIDCut,
	    gammaCuts.passGammaPCALFiducialCut,
	]
    }

    def setCutParameterFromMagField( field_config ){
		gammaCuts.setGammaCutParameters(field_config)
    }

    def setGammaCutStrictness(gam_cut_strictness){
		gammaCuts.setGammaCutStrictness(gam_cut_strictness)
    }


    def getGoodGammaWithCuts(event){
	//return map - key is index of track in REC::Particle and value is list of booleans for the cuts
	def gam_cut_result = (0..<event.npart).findAll{event.charge[it]==0}.collect{ ii -> [ii, gammaCutStrategies.collect{ gam_test -> gam_test(event,ii) } ] }.collectEntries()	  		
	return gam_cut_result
    }

    def getGoodGamma(event){
	//return a list of REC::Particle indices for tracks passing all gamma cuts
	def gam_cut_result = (0..<event.npart).findAll{event.charge[it]==0}.collect{ ii -> [ii, gammaCutStrategies.collect{ gam_test -> gam_test(event,ii) } ] }.collectEntries()
	return gam_cut_result.findResults{gam_indx, cut_result -> !cut_result.contains(false) ? gam_indx : null}
    }

}
