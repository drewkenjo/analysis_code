import java.io.*;
import java.util.*;
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.graphics.EmbeddedCanvas;
import pid.sangbaek.electron

def run = args[0].toInteger()

TDirectory out = new TDirectory()
out.mkdir('/elec/')
out.cd('/elec/')

def EB = 10.6f
if(run>6607) EB=10.2f

def binnum = args[1].toInteger()

HistoDef(EB, binnum)

filenum=-1

for (arg in args){
	filenum++
	if (filenum<2) continue
	HipoDataSource reader = new HipoDataSource();
	reader.open(arg);
	while( reader.hasEvent()){
		DataEvent event = reader.getNextEvent();
		processEvent(event)
	}
}

(0..<6).each{
	out.addDataSet(H_elec_vz[it])
	out.addDataSet(H_elec_nphe[it])
	out.addDataSet(H_elec_Sampl[it])
	out.addDataSet(H_elec_vz_mom[it])
	out.addDataSet(H_elec_nphe_mom[it])
	out.addDataSet(H_elec_Sampl_mom[it])
	out.addDataSet(H_elec_PCALECAL[it])
	out.addDataSet(H_neg_vz[it])
	out.addDataSet(H_neg_nphe[it])
	out.addDataSet(H_neg_Sampl[it])
	out.addDataSet(H_neg_vz_mom[it])
	out.addDataSet(H_neg_nphe_mom[it])
	out.addDataSet(H_neg_Sampl_mom[it])
	out.addDataSet(H_neg_PCALECAL[it])
}

out.writeFile('elec_'+run+'_bin_'+binnum+'.hipo')

public void processEvent(DataEvent event) {
	if(!event.hasBank("REC::Particle")) return 
	if(!event.hasBank("REC::Calorimeter")) return 
    def evc = event.getBank("REC::Calorimeter")
    def secs = [evc.getShort('pindex')*.toInteger(), evc.getByte('sector')].transpose().collectEntries()
    def evp = event.getBank("REC::Particle")
    def eleind = electron.find_byBANK(evp)
    evp.getInt("pid").eachWithIndex{pid, ind ->
    	if(secs[ind]==null) return
	    int charge = evp.getInt("charge",ind)
	    vz = evp.getFloat("vz",ind)
	    mom = (float) ['x','y','z'].collect{evp.getFloat("p"+it,ind)}.sum()
	    e_ecal_E = 0
	    e_etot_E = 0
	    e_pcal_E = 0
	    evc.getInt("pindex").eachWithIndex{ pindex, ind_c ->
	    	if (pindex==ind){
	    		det = evc.getInt("layer", ind_c);
	    		e_ecal_E+=evc.getFloat("energy",ind_c)
	    		if (det==1) e_pcal_E+=evc.getFloat("energy",ind_c)
	    		else if (det==4 || det==7) e_etot_E+=evc.getFloat("energy",ind_c)
	    	}
	    }
   	    sampl_frac = e_ecal_E/mom

   	    // Negative
   	    if (charge<0){
	    	H_neg_vz[secs[ind]-1].fill(vz) //vz
	    	if (e_ecal_E>0) {
	    		H_neg_Sampl[secs[ind]-1].fill(sampl_frac) // sampling Fraction
		    	H_neg_Sampl_mom[secs[ind]-1].fill(mom,sampl_frac)
		    }
	    	if (e_pcal_E>0 && e_etot_E>0) H_neg_PCALECAL[secs[ind]-1].fill(e_pcal_E,e_etot_E);
	    	H_neg_vz_mom[secs[ind]-1].fill(mom,vz)
	    	if(!event.hasBank("REC::Cherenkov")) return
	    	def evh = event.getBank("REC::Cherenkov")
	    	evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
	    		if(evh.getInt("detector",ind_h)!=15) return
    			if(pindex==ind) {
    				H_neg_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
    				H_neg_nphe_mom[secs[ind]-1].fill(mom,evh.getFloat("nphe",ind_h))
    			}
	    	}
    	}

    	// Electron
    	if (ind==eleind){
    		H_elec_vz[secs[ind]-1].fill(vz)
	    	H_elec_vz_mom[secs[ind]-1].fill(mom,vz)
	    	if (e_ecal_E>0) {
	    		H_elec_Sampl[secs[ind]-1].fill(sampl_frac) // sampling Fraction
		    	H_elec_Sampl_mom[secs[ind]-1].fill(mom,sampl_frac)
		    }
	    	if (e_pcal_E>0 && e_etot_E>0) H_elec_PCALECAL[secs[ind]-1].fill(e_pcal_E,e_etot_E);
	    	if(!event.hasBank("REC::Cherenkov")) return
	    	def evh = event.getBank("REC::Cherenkov")
	    	evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
	    		if(evh.getInt("detector",ind_h)!=15) return
    			if(pindex==ind) {
					H_elec_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
					H_elec_nphe_mom[secs[ind]-1].fill(mom,evh.getFloat("nphe",ind_h))
    			}
	    	}
    	}
    }
}

public void HistoDef(float eb, int binnum){
	H_elec_vz =(0..<6).collect{new H1F("H_elec_vz_S"+(it+1), "H_elec_vz_S"+(it+1),binnum,-15,15);}
	H_elec_nphe =(0..<6).collect{new H1F("H_elec_nphe_S"+(it+1), "H_elec_nphe_S"+(it+1),binnum,0,50);}
	H_elec_Sampl =(0..<6).collect{new H1F("H_elec_Sampl_S"+(it+1), "H_elec_Sampl_S"+(it+1),binnum,0,1);}
	H_elec_vz_mom =(0..<6).collect{new H2F("H_elec_vz_mom_S"+(it+1), "H_elec_vz_mom_S"+(it+1),binnum,0,eb,binnum,-15,15);}
	H_elec_nphe_mom =(0..<6).collect{new H2F("H_elec_nphe_mom_S"+(it+1), "H_elec_nphe_mom_S"+(it+1),binnum,0,eb,binnum,0,50);}
	H_elec_Sampl_mom =(0..<6).collect{new H2F("H_elec_Sampl_mom_S"+(it+1), "H_elec_Sampl_mom_S"+(it+1),binnum,0,eb,binnum,0,1);}
	H_elec_PCALECAL = (0..6).collect{new H2F("H_elec_PCALECAL_S"+(it+1),"H_elec_PCALECAL_S"+(it+1),100,0,1.5,100,0,1.5);}

	H_neg_vz =(0..<6).collect{new H1F("H_neg_vz_S"+(it+1), "H_neg_vz_S"+(it+1),binnum,-15,15);}
	H_neg_nphe =(0..<6).collect{new H1F("H_neg_nphe_S"+(it+1), "H_neg_nphe_S"+(it+1),binnum,0,50);}
	H_neg_Sampl =(0..<6).collect{new H1F("H_neg_Sampl_S"+(it+1), "H_neg_Sampl_S"+(it+1),binnum,0,1);}
	H_neg_vz_mom =(0..<6).collect{new H2F("H_neg_vz_mom_S"+(it+1), "H_neg_vz_mom_S"+(it+1),binnum,0,eb,binnum,-15,15);}
	H_neg_nphe_mom =(0..<6).collect{new H2F("H_neg_nphe_mom_S"+(it+1), "H_neg_nphe_mom_S"+(it+1),binnum,0,eb,binnum,0,50);}
	H_neg_Sampl_mom =(0..<6).collect{new H2F("H_neg_Sampl_mom_S"+(it+1), "H_neg_Sampl_mom_S"+(it+1),binnum,0,eb,binnum,0,1);}
	H_neg_PCALECAL = (0..6).collect{new H2F("H_neg_PCALECAL_S"+(it+1), "H_neg_PCALECAL_S"+(it+1),100,0,1.5,100,0,1.5);}
}
