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
import electron

def run = args[0].toInteger()

TDirectory out = new TDirectory()
out.mkdir('/elec/')
out.cd('/elec/')

def EB = 10.6f
if(run>6607) EB=10.2f

def binnum = args[1].toInteger()

HistoDef(EB, binnum)

int e_part_ind, e_sect, e_nphe
float e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_ecal_E, e_Sampl_frac
LorentzVector Ve = new LorentzVector()

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
	out.addDataSet(H_elec_HTCC_nphe[it])
	out.addDataSet(H_elec_EC_Sampl[it])
	out.addDataSet(H_elec_mom_nphe[it])
	out.addDataSet(H_elec_mom_Sampl[it])
	out.addDataSet(H_neg_vz[it])
	out.addDataSet(H_neg_HTCC_nphe[it])
	out.addDataSet(H_neg_EC_Sampl[it])
	out.addDataSet(H_neg_mom_nphe[it])
	out.addDataSet(H_neg_mom_Sampl[it])
}

out.writeFile('electron_pID_'+run+'_bin_'+binnum+'.hipo')

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
	    energy = 0
	    evc.getInt("pindex").eachWithIndex{ pindex, ind_c ->
	    	if (pindex==ind) energy+=evc.getFloat("energy",ind_c)
	    }
   	    sampl_frac = energy/mom

   	    // Negative
   	    if (charge<0){
	    	H_neg_vz[secs[ind]-1].fill(vz) //vz
	    	H_neg_EC_Sampl[secs[ind]-1].fill(sampl_frac) // sampling Fraction
	    	H_neg_mom_Sampl[secs[ind]-1].fill(sampl_frac,mom)
	    	if(!event.hasBank("REC::Cherenkov")) return
	    	def evh = event.getBank("REC::Cherenkov")
	    	evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
	    		if(evh.getInt("detector",ind_h)!=15) return
    			if(pindex==ind) {
    				H_neg_HTCC_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
    				H_neg_mom_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h),mom)
    			}
	    	}
    	}

    	// Electron
    	if (ind==eleind){
    		H_elec_vz[secs[ind]-1].fill(vz)
			H_elec_EC_Sampl[secs[ind]-1].fill(sampl_frac) //electron
	    	H_elec_mom_Sampl[secs[ind]-1].fill(sampl_frac,mom)
	    	if(!event.hasBank("REC::Cherenkov")) return
	    	def evh = event.getBank("REC::Cherenkov")
	    	evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
	    		if(evh.getInt("detector",ind_h)!=15) return
    			if(pindex==ind) {
					H_elec_HTCC_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
					H_elec_mom_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h),mom)
    			}
	    	}
    	}
    }
}

public void HistoDef(float eb, int binnum){
	H_elec_vz =(0..<6).collect{new H1F("H_elec_vz_S"+(it+1), "H_elec_vz_S"+(it+1),binnum,-25,25);}
	H_elec_HTCC_nphe =(0..<6).collect{new H1F("H_elec_HTCC_nphe_S"+(it+1), "H_elec_HTCC_nphe_S"+(it+1),binnum,0,100);}
	H_elec_EC_Sampl =(0..<6).collect{new H1F("H_elec_EC_Sampl_S"+(it+1), "H_elec_EC_Sampl_S"+(it+1),binnum,0,1);}
	H_elec_mom_nphe =(0..<6).collect{new H2F("H_elec_mom_nphe_S"+(it+1), "H_elec_mom_nphe_S"+(it+1),binnum,0,100,binnum,0,eb);}
	H_elec_mom_Sampl =(0..<6).collect{new H2F("H_elec_mom_Sampl_S"+(it+1), "H_elec_mom_Sampl_S"+(it+1),binnum,0,1,binnum,0,eb);}
	H_neg_vz =(0..<6).collect{new H1F("H_neg_vz_S"+(it+1), "H_neg_vz_S"+(it+1),binnum,-25,25);}
	H_neg_HTCC_nphe =(0..<6).collect{new H1F("H_neg_HTCC_nphe_S"+(it+1), "H_neg_HTCC_nphe_S"+(it+1),binnum,0,100);}
	H_neg_EC_Sampl =(0..<6).collect{new H1F("H_neg_EC_Sampl_S"+(it+1), "H_neg_EC_Sampl_S"+(it+1),binnum,0,1);}
	H_neg_mom_nphe =(0..<6).collect{new H2F("H_neg_mom_nphe_S"+(it+1), "H_neg_mom_nphe_S"+(it+1),binnum,0,100,binnum,0,eb);}
	H_neg_mom_Sampl =(0..<6).collect{new H2F("H_neg_mom_Sampl_S"+(it+1), "H_neg_mom_Sampl_S"+(it+1),binnum,0,1,binnum,0,eb);}
}