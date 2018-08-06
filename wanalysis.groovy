import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;


double en = Double.parseDouble(args[1]);
double enmax = en+0.1; //GeV
double thetamax = 40;  //degrees
double phimax = 3.2;   //radians
double vzmax = 50;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

HipoDataSource reader = new HipoDataSource();

H1F momentum = new H1F("momentum", "momentum", 500, 0, 10);
momentum.setTitleX("momentum");

H1F W_hist = new H1F("W", "W", 500, 0, wmax);
W_hist.setTitleX("W");

H1F Q2_hist = new H1F("Q2", "Q2", 500, 0, enmax);
Q2_hist.setTitleX("Q2");

H2F W_vs_Q2 = new H2F("W_vs_Q2", "W_vs_Q2", 500, 0.0, enmax, 500, 0.0, wmax);
W_vs_Q2.setTitleX("Q2");
W_vs_Q2.setTitleY("W");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E_vs_Theta", 500, 5, thetamax, 500, 0, enmax);
E_vs_Theta.setTitleX("Theta");
E_vs_Theta.setTitleY("E'");

H2F z_vs_Theta = new H2F("z_vs_Theta", "z_vs_Theta", 500, 5, thetamax, 500, -vzmax, vzmax);
z_vs_Theta.setTitleX("Theta");
z_vs_Theta.setTitleY("z vertex");

H2F Phi_vs_Theta = new H2F("Phi_vs_Theta", "Phi_vs_Theta", 500, 5, thetamax, 500, -phimax, phimax);
Phi_vs_Theta.setTitleX("Theta");
Phi_vs_Theta.setTitleY("Phi");

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 500, 0, wmax, 500, -phimax, phimax);
Phi_vs_W.setTitleX("W");
Phi_vs_W.setTitleY("Phi");
	
double e_mass = 0.000511;
double p_mass = 0.93827203;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);
reader.open(args[0]);

double emax = 0;
phimax = 0;
thetamax = 0;
vzmax = 0;
int counter = 0;

while (reader.hasEvent()) {
	DataEvent event = reader.getNextEvent();
      	if (event.hasBank("REC::Particle")) {
      	 	DataBank bank_rec = event.getBank("REC::Particle");
	 		//counter++;
	 		//if(counter > 100){break;}
		for (int k = 0; k < bank_rec.rows(); k++) {
			int pid = bank_rec.getInt("pid", k);
			float px = bank_rec.getFloat("px", k);
			float py = bank_rec.getFloat("py", k);								
			float pz = bank_rec.getFloat("pz", k);
			float beta = bank_rec.getFloat("beta", k);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			double phi = Math.atan2((double) py,(double) px);
			double theta = Math.acos((double) pz/(double) mom);
			theta *= 180/Math.PI;
			float vz = bank_rec.getFloat("vz", k);	     
			
			if (pid != 11) continue;
			
			momentum.fill(mom);
			
			Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
			LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
			e_vec_prime.setVectM(e_vec_3, e_mass);
			
			if(e_vec_prime.e() < 0.1 * en){continue;}
		
			LorentzVector q_vec = new LorentzVector(); //4 vector q
			q_vec.copy(e_vec); //e - e'
			q_vec.sub(e_vec_prime);
			double Q2 = -q_vec.mass2(); //-q^2
			Q2_hist.fill(Q2);
			
			LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
			w_vec.copy(p_vec); //p-q
			w_vec.add(q_vec);
			double W = w_vec.mass(); 
			W_hist.fill(W);
			W_vs_Q2.fill(Q2,W);
			Phi_vs_W.fill(W,phi);
			
			if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} //calculate max values of each param
			if(theta > thetamax){thetamax = theta;}
			if(phi > phimax){phimax = phi;}
			if(vz > vzmax){vzmax = vz;}
			if(theta > 5)
			{
				E_vs_Theta.fill(theta,e_vec_prime.e());
				z_vs_Theta.fill(theta,vz);
				Phi_vs_Theta.fill(theta,phi);
			}
		}
	}
}

System.out.println(emax + " " + thetamax + " " + phimax + " " + vzmax);

TCanvas can = new TCanvas("can", 800, 600);
can.draw(W_vs_Q2);
can.save("W_vs_Q2.png");

TCanvas can1 = new TCanvas("can", 800, 600);
can1.draw(momentum);
can1.save("mom.png");

TCanvas can2 = new TCanvas("can", 800, 600);
can2.draw(Q2_hist);
can2.save("Q2.png");

TCanvas can3 = new TCanvas("can", 800, 600);
can3.draw(W_hist);
can3.save("W.png");

TCanvas can4 = new TCanvas("can", 800, 600);
can4.draw(E_vs_Theta);
can4.save("EvsTheta.png");

TCanvas can5 = new TCanvas("can", 800, 600);
can5.draw(z_vs_Theta);
can5.save("ZvsTheta.png");

TCanvas can6 = new TCanvas("can", 800, 600);
can6.draw(Phi_vs_Theta);
can6.save("PhivsTheta.png");

TCanvas can7 = new TCanvas("can", 800, 600);
can7.draw(Phi_vs_W);
can7.save("PhivsW.png");
