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
double phimax = 180;   //degrees
double vzmax = 50;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

System.out.println("Energy is " + en + " W max is " + wmax);

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

H2F Cal_y_vs_x_precut = new H2F("Cal_y_vs_x_precut", "Cal_y_vs_x_precut", 500, -450,450, 500, -450,450);
Cal_y_vs_x_precut.setTitleX("X (cm)");
Cal_y_vs_x_precut.setTitleY("Y (cm)");

H1F Cal_lu = new H1F("Cal_lu", "Cal_lu", 500, 0, 450); 
H1F Cal_lv = new H1F("Cal_lv", "Cal_lv", 500, 0, 450); 
H1F Cal_lw = new H1F("Cal_lw", "Cal_lw", 500, 0, 450); 

H2F Cal_y_vs_x = new H2F("Cal_y_vs_x", "Cal_y_vs_x", 500, -450,450, 500, -450, 450);
Cal_y_vs_x.setTitleX("X (cm)");
Cal_y_vs_x.setTitleY("Y (cm)");

H2F DCXvsHTCCX = new H2F("DCXvsHTCCX", "DCXvsHTCCX", 500, -450,450,500,-50,50);
DCXvsHTCCX.setTitleX("DC X");
DCXvsHTCCX.setTitleY("delta X");

H2F DCYvsHTCCY = new H2F("DCYvsHTCCY", "DCYvsHTCCY", 500, -450,450,500,-50,50);
DCYvsHTCCY.setTitleX("DC Y");
DCYvsHTCCY.setTitleY("delta Y");


H2F DCPhivsHTCCPhi = new H2F("DCPhivsHTCCPhi", "DCPhivsHTCCPhi", 500, -180,180,500,-10,10);
DCPhivsHTCCPhi.setTitleX("DC Phi");
DCPhivsHTCCPhi.setTitleY("delta Phi");

H2F PhiPvsPhiE = new H2F("PhiPvsPhiE","PhiPvsPhiE",500, -180, 180, 500, -180, 180);
PhiPvsPhiE.setTitleX("Phi Electron");
PhiPvsPhiE.setTitleY("Phi Proton");

H2F EnPvsEnE = new H2F("EnPvsEnE","EnPvsEnE",500, 0 , enmax, 500, 0 , enmax);
EnPvsEnE.setTitleX("Energy Electron");
EnPvsEnE.setTitleY("Energy Proton");

H2F ZPvsZE = new H2F("ZPvsZE","ZPvsZE",500, 0 , enmax, 500, 0 , enmax);
ZPvsZE.setTitleX("vertex Electron");
ZPvsZE.setTitleY("vertex Proton");

H1F PhiResPE = new H1F("PhiResPE","PhiResPE", 500, 0, 360);
PhiResPE.setTitleX("Phi p - Phi e");

H1F EResPE = new H1F("EResPE","EResPE", 500, -10, 10);
EResPE.setTitleX("E p - E e");

H1F VResPE = new H1F("VResPE","VResPE", 500, -10, 10);
VResPE.setTitleX("Vz P - Vz E");



HashMap<Integer,H1F> histmap = new HashMap<Integer,H1F>();
for(int i = 5; i <= 20; i++){histmap.put(i,new H1F("Phi vs Theta " + i, 500,-phimax,phimax));}
	
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
byte sector = 0;
int cal_row = 0;
int dc_row = 0;
int htcc_row = 0;
float x_dc = 1000.0;
float y_dc = 1000.0;
float z_dc = 1000.0;
float x_htcc = 1000.0;
float y_htcc = 1000.0;
float z_htcc = 1000.0;
double phi_dc = 1000.0;
double phi_htcc = 1000.0;
float e_px = 0;
float e_py = 0;
float e_pz = 0; 
float e_beta = 0;
float e_vz = 0;
float e_phi = 0;
float e_mom = 0;
float e_theta = 0;

float p_px = 0;
float p_py = 0;
float p_pz = 0; 
float p_beta = 0;
float p_vz = 0;
float p_phi = 0;
float p_theta = 0; 
float p_mom = 0;
 
boolean found_electron = false;
boolean found_proton = false; 

Vector3 e_vec_3 = new Vector3();
Vector3 p_vec_3 = new Vector3();
LorentzVector e_vec_prime = new LorentzVector();
LorentzVector p_vec_prime = new LorentzVector();

while (reader.hasEvent()) {
	DataEvent event = reader.getNextEvent();
	if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Traj")) {
		DataBank bank_rec = event.getBank("REC::Particle");
		DataBank bank_cal = event.getBank("REC::Calorimeter");
		DataBank bank_traj = event.getBank("REC::Traj");
			counter++;
			//if(counter > 400){break;}
		for (int k = 0; k < bank_rec.rows(); k++) {
			int pid = bank_rec.getInt("pid", k);
			byte q = bank_rec.getByte("charge", k);
			float px = bank_rec.getFloat("px", k);
			float py = bank_rec.getFloat("py", k);								
			float pz = bank_rec.getFloat("pz", k);
			float beta = bank_rec.getFloat("beta", k);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			double phi = Math.atan2((double) py,(double) px);
			double theta = Math.acos((double) pz/(double) mom);
			theta *= 180/Math.PI;
			phi *= 180/Math.PI;
			float vz = bank_rec.getFloat("vz", k);	     

			//if (pid != 11) continue;
			//if(q != -1) continue;

			if(pid == 11){
			        found_electron = true;
				e_px = px; 
				e_py = py; 
				e_pz = pz; 
				e_beta = beta;
				e_mom = mom; 
				e_phi = phi; 
				e_theta = theta; 
				e_vz = vz; 
				
				e_vec_3 = new Vector3(px, py, pz); //3 vector e'
				e_vec_prime = new LorentzVector(); //4 vector e'
				e_vec_prime.setVectM(e_vec_3, e_mass);
				if(e_vec_prime.e() < 0.1 * en){continue;} //cut below 10% beam
				if(theta < 5 || theta > 40){continue;} //cut outside of 5 and 40 degrees for FD

				cal_row = cal_cut_row(event, k);
				//System.out.println(j + " " + bank_cal.rows());
				if(cal_row != -1){
					sector = bank_cal.getByte("sector",cal_row);
					float x_cal = bank_cal.getFloat("x",cal_row);
					float y_cal = bank_cal.getFloat("y",cal_row);
					float lu = bank_cal.getFloat("lu",cal_row);
					float lv = bank_cal.getFloat("lv",cal_row);
					float lw = bank_cal.getFloat("lw",cal_row);
					Cal_y_vs_x_precut.fill(x_cal,y_cal);
					if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
						Cal_lu.fill(lu);
						Cal_lv.fill(lv);
						Cal_lw.fill(lw);
						Cal_y_vs_x.fill(x_cal,y_cal);
					}
				}

				dc_row = dc_cut_row(event, k);
				if(dc_row != -1){
					x_dc = bank_traj.getFloat("x",dc_row);
					y_dc = bank_traj.getFloat("y",dc_row);
					z_dc = bank_traj.getFloat("z",dc_row);
					double pos_dc = Math.sqrt(x_dc*x_dc + y_dc*y_dc + z_dc*z_dc);
					double theta_dc = Math.acos((double) z_dc/ pos_dc);
					phi_dc = Math.atan2((double) y_dc,(double) x_dc);
					theta_dc *= 180/Math.PI;
					phi_dc *= 180/Math.PI;
					//if(dc_cut(x_dc,y_dc,sector)){
						//System.out.println(theta_dc);
						if(histmap.containsKey((int) Math.floor(theta_dc))){
							histmap.get((int) Math.floor(theta_dc)).fill(phi_dc);
						}
					//}
				}
				htcc_row = htcc_cut_row(event,k);
				if(htcc_row != -1){
					x_htcc = bank_traj.getFloat("x",htcc_row);
					y_htcc = bank_traj.getFloat("y",htcc_row);
					z_htcc = bank_traj.getFloat("z",htcc_row);
					double pos_htcc = Math.sqrt(x_htcc*x_htcc + y_htcc*y_htcc + z_htcc*z_htcc);
					double theta_htcc = Math.acos((double) z_htcc/ pos_htcc);
					phi_htcc = Math.atan2((double) y_htcc,(double) x_htcc);
					theta_htcc *= 180/Math.PI;
					phi_htcc *= 180/Math.PI;
				}
				if(dc_row != -1 && htcc_row != -1)
				{
					//if(counter%1000 == 1) System.out.println(dc_row + " " + htcc_row);

					if(x_htcc != 1000.0 && x_dc != 1000.0)
					{	  
						  DCXvsHTCCX.fill(x_dc,x_dc-x_htcc);
					}
					if(y_htcc != 1000.0 && y_dc != 1000.0)
					{	  
						  DCYvsHTCCY.fill(y_dc,y_dc-y_htcc);
					}
					if(phi_htcc != 1000.0 && phi_dc != 1000.0)
					{					  
						  DCPhivsHTCCPhi.fill(phi_dc ,phi_dc-phi_htcc);
					}	  
				}

				momentum.fill(mom);
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

				E_vs_Theta.fill(theta,e_vec_prime.e());
				z_vs_Theta.fill(theta,vz);
				Phi_vs_Theta.fill(theta,phi);
			}
			if(pid == 2212){
			       	found_proton = true;
				p_px = px; 
				p_py = py; 
				p_pz = pz; 
				p_beta = beta;
				p_mom = mom; 
				p_phi = phi; 
				p_theta = theta; 
				p_vz = vz; 
				
				p_vec_3 = new Vector3(px, py, pz); //3 vector e'
				p_vec_prime = new LorentzVector(); //4 vector e'
				p_vec_prime.setVectM(p_vec_3, p_mass);
			}



		}
		if(found_electron && found_proton){
			PhiPvsPhiE.fill(e_phi,p_phi);
			EnPvsEnE.fill(e_vec_prime.e(),p_vec_prime.e());
			ZPvsZE.fill(e_vz,p_vz);
			PhiResPE.fill(Math.abs(e_phi-p_phi));
			EResPE.fill(e_vec_prime.e()-p_vec_prime.e());
			VResPE.fill(e_vz-p_vz);
		}
		found_electron = false; 
		found_proton = false;
		
	}
	/*if(event.hasBank("RECHB::Calorimeter")){
		DataBank bank_cal = event.getBank("RECHB::Calorimeter");
		for(int j = 0; j < bank_cal.rows(); j++){
			float x = bank_cal.getFloat("x",j);
			float y = bank_cal.getFloat("y",j);
			float lu = bank_cal.getFloat("lu",j);
			float lv = bank_cal.getFloat("lv",j);
			float lw = bank_cal.getFloat("lw",j);
			Cal_y_vs_x_precut.fill(x,y);
			if(lu > 350 || lu < 60 || lv > 370 || lw > 390){continue;}
			Cal_lu.fill(lu);
			Cal_lv.fill(lv);
			Cal_lw.fill(lw);
			Cal_y_vs_x.fill(x,y);
		}	
	}*/
}


boolean dc_cut(float X, float Y, int S)
{ 
	boolean result= false;
	if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
	&& (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
  
	return result;
}

int cal_cut_row(DataEvent event, int row){
	DataBank bank_cal = event.getBank("REC::Calorimeter");
	int row_index = 0;
	int cal_row_match = -1;
	for(int j = 0; j < bank_cal.rows(); j++){
		row_index = bank_cal.getInt("pindex",j);
		if(row_index == row){
			cal_row_match = j;
			break;
		}
	}
	return cal_row_match;
}

int dc_cut_row(DataEvent event, int row){
	DataBank bank_traj = event.getBank("REC::Traj");
	int row_index = 0;
	int det_id = 0;
	int cal_row_match = -1;
	float path_length = 0;
	for(int j = 0; j < bank_traj.rows(); j++){
		row_index = bank_traj.getInt("pindex",j);
		det_id = bank_traj.getInt("detId",j);
		path_length = bank_traj.getFloat("pathlength",j);
		if(row_index == row && path_length > 250 && path_length < 550){
			cal_row_match = j;
			break;
		}
	}
	return cal_row_match;
}
int htcc_cut_row(DataEvent event, int row){
	DataBank bank_traj = event.getBank("REC::Traj");
	int row_index = 0;
	int det_id = 0;
	int htcc_row_match = -1;
	float path_length = 0;
	for(int j = 0; j < bank_traj.rows(); j++){
		row_index = bank_traj.getInt("pindex",j);
		det_id = bank_traj.getInt("detId",j);
		//path_length = bank_traj.getFloat("pathlength",j);
		if(row_index == row && det_id == 0){
			htcc_row_match = j;
			break;
		}
	}
	return htcc_row_match;
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

TCanvas can8 = new TCanvas("can", 800,600);
can8.draw(Cal_y_vs_x_precut);
can8.save("Calyvxprecut.png");
	   
TCanvas can9 = new TCanvas("can", 800,600);
can9.draw(Cal_lu);
can9.save("Callu.png");
	   
TCanvas can10 = new TCanvas("can", 800,600);
can10.draw(Cal_lv);
can10.save("Callv.png");
	   
TCanvas can11 = new TCanvas("can", 800,600);
can11.draw(Cal_lw);
can11.save("Callw.png");

TCanvas can12 = new TCanvas("can", 800,600);
can12.draw(Cal_y_vs_x);
can12.save("Calyvx.png");
	   
TCanvas can13 = new TCanvas("can", 800,600);
can13.draw(DCXvsHTCCX);
can13.save("DCXvsHTCCX.png");

TCanvas can14 = new TCanvas("can", 800,600);
can14.draw(DCPhivsHTCCPhi);
can14.save("DCPhivsHTCCPhi.png");

TCanvas can15 = new TCanvas("can", 800,600);
can15.draw(DCYvsHTCCY);
can15.save("DCYvsHTCCY.png");

TCanvas can16 = new TCanvas("can", 800,600);
can16.draw(PhiPvsPhiE);
can16.save("PhiPvsPhiE.png");

TCanvas can17 = new TCanvas("can", 800,600);
can17.draw(EnPvsEnE);
can17.save("EnPvsEnE.png");

TCanvas can18 = new TCanvas("can", 800,600);
can18.draw(ZPvsZE);
can18.save("ZPvsZE.png");

TCanvas can19 = new TCanvas("can", 800,600);
can19.draw(PhiResPE);
can19.save("PhiResPE.png");

TCanvas can20 = new TCanvas("can", 800,600);
can20.draw(EResPE);
can20.save("EResPE.png");

TCanvas can21 = new TCanvas("can", 800,600);
can21.draw(VResPE);
can21.save("VResPE.png");

HashMap<Integer,TCanvas> canvasmap = new HashMap<Integer,TCanvas>();
for(int i : histmap.keySet()){
	canvasmap.put(i, new TCanvas("can", 800,600));
	canvasmap.get(i).draw(histmap.get(i));
	canvasmap.get(i).save("phivstheta" + i + ".png");
}
	   
System.out.println("Done!");
