#!/opt/coatjava/6.3.1/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.jroot.ROOTFile
//groovy n java
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
//personal 
import examples.epkpkm
import mon.electronPlotter
import event.EventConverter
import event.Event

//dmriser
LorentzVector.metaClass.static.withPID = {pid,x,y,z->
    double m = PDGDatabase.getParticleMass(pid)
    double e = Math.sqrt(x*x+y*y+z*z+m*m)
    return new LorentzVector(x,y,z,e)
}

LorentzVector.metaClass.plus = {
    LorentzVector a = new LorentzVector(delegate)
    a.add(it)
    return a
}

LorentzVector.metaClass.minus = {
    LorentzVector a = new LorentzVector(delegate)
    a.sub(it)
    return a
}

def outname = args[0].split('/')[-1]
def out = new TDirectory()

def processors = [new electronPlotter()]
processors.each{out.mkdir("/"+it.getClass().getSimpleName())}

def timer = new Timer()

timer.schedule({
    processors.each{
	out.cd("/"+it.getClass().getSimpleName())
	it.histos.each{out.addDataSet(it.value)}
    }
    out.writeFile("test_$outname")
} as TimerTask, 30000,60000)


GParsPool.withPool 16, {
    args.eachParallel{fname->
	def reader = new HipoDataSource()
	reader.open(fname)
	
	while(reader.hasEvent()) {
	    def data_event = reader.getNextEvent()
	    
	    def event = EventConverter.convert(data_event)
	    
	    processors.each{it.processEvent(event)}

	}

	reader.close()
    }
}


timer.cancel()

//finally write output to a ROOT file
def outroot = new ROOTFile("electronPID_$outname"+".root")
processors.each{channel -> channel.histos.each{outroot.writeDataSet(it.value)}}
outroot.close()
