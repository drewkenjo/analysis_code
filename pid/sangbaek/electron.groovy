import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

class electron{
	static def find_byPID = { pid ->
    pid==11
  }

  static def find_bySTATUS = { status ->
    status<0
  }

  static def find_byCHARGE = {charge ->
    charge<0
  }

  static def find_byNPhe = {nphe ->
    nphe>2
  }

  static def find_bySampl = {sampl ->
    sampl >0.18
  }

  static def find_byVZ = {vz ->
    Math.abs(vz) < 1
  }
  static def find_byBANK = {pbank ->
     return (0..pbank.rows())
      .find{ ind->
        def pid = pbank.getInt('pid',ind)
        def status = pbank.getShort('status',ind)
        def charge = pbank.getByte('charge',ind)
        this.find_byPID(pid) && this.find_bySTATUS(status) && this.find_byCHARGE(charge)
      }
  }

  static def find_byEVENT = { event ->

    def pbank = event.getBank("REC::Particle")
    return this.find_byBANK(pbank)
  }
}