package utils
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase


class KinTool{

    //move somplace else

    static def calcQ2(LorentzVector beam, LorentzVector measured_el){
        LorentzVector VGS = beam - measured_el
	return -VGS.mass2()
    }

    static def calcXb(LorentzVector beam, LorentzVector measured_el){
        LorentzVector VGS = beam - measured_el
	return calcQ2(beam,measured_el)/(2.0*PDGDatabase.getParticleMass(2212)*VGS.e())
    }

    static def calcT(LorentzVector measured_prot){
	return 2.0*PDGDatabase.getParticleMass(2212)*(measured_prot.e()-PDGDatabase.getParticleMass(2212))
    }

    static def calcT2(LorentzVector measured_el, LorentzVector measured_gam){
    return 2.0*PDGDatabase.getParticleMass(2212)*(measured_prot.e()-PDGDatabase.getParticleMass(2212))
    }

    //used for DVCS
    static def calcT2(LorentzVector beam, LorentzVector measured_el, LorentzVector measured_gam){
        def Q2 = calcQ2(beam, measured_el)
        def nu = calcNu(beam, measured_el)
        LorentzVector VGS = beam - measured_el
        def M = PDGDatabase.getParticleMass(2212)
        def costheta = VGS.vect().dot(measured_gam.vect())/VGS.vect().mag()/measured_gam.vect().mag()
        return (M*Q2+2*M*nu*(nu-Math.sqrt(nu*nu+Q2)*costheta))/(M+nu-Math.sqrt(nu*nu+Q2)*costheta)
    }


    static def calcNu(LorentzVector beam, LorentzVector measured_el){
	return calcQ2(beam,measured_el)/(2*PDGDatabase.getParticleMass(2212)*calcXb(beam,measured_el))
    }
    
    static def calcY(LorentzVector beam, LorentzVector measured_el){
	return calcNu(beam,measured_el)/beam.e()
    }

    static def calcPhiTrento(LorentzVector beam, LorentzVector measured_el,LorentzVector measured_prot){
    def v3l = beam.vect().cross(measured_el.vect())
    LorentzVector VGS = beam - measured_el
    def v3h = measured_prot.vect().cross(VGS.vect())
    def trento = Vangle(v3l, v3h)

    if(v3l.dot(measured_prot.vect())>0){trento=360-trento}
    return trento
    }

    // used for DVCS.
    static def calcPhiTrento2(LorentzVector beam, LorentzVector measured_el,LorentzVector measured_gam){
    LorentzVector VGS = beam - measured_el
    def v1 = VGS.vect().cross(measured_el.vect())
    def v2 = VGS.vect().cross(measured_gam.vect())
    def trento = KinTool.Vangle(v1,v2)

    if (VGS.vect().dot(v1.cross(v2))<0) trento = 360 - trento
    return trento
    }


    static def delta_meas_energy( Double beam,  LorentzVector measured_el ){
	def calc_e = beam/(1+ (beam/PDGDatabase.getParticleMass(2212))*(1-Math.cos(measured_el.theta())) );
	return (calc_e  - measured_el.e() )
    }

    static def delta_meas_theta( Double beam,  LorentzVector measured_el ){
	def calc_theta = Math.toDegrees(Math.acos( 1 + (PDGDatabase.getParticleMass(2212)/beam)*( 1 - beam/measured_el.e()) ))
	return Math.toDegrees((calc_theta  - measured_el.theta() ))
    }

    static def Vangle = {v1, v2 -> 
        if( v1.mag() * v2.mag() !=0 && v1.dot(v2)<v1.mag()*v2.mag() ) return Math.toDegrees( Math.acos(v1.dot(v2)/(v1.mag()*v2.mag()) ) ); 
        else return 0
    }
    
}


