/**
 * 
 */
package src;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import src.Constant.VALUE_TYPE;

/**
 * @author dream
 *
 */
public class ClusterValidityIndices {
	private ArrayList<Neuron> listNeuron = null;
	private ExprParm para;
	private ArrayList<InputData> ListTrnInput;
	
	public ClusterValidityIndices(ExprParm para){
		this.para = para;
		ListTrnInput = para.getListTrnInput();
	}
	
	/**
	 * 
	 * @param vtCluster
	 * @return
	 */
	/*
	public Vector<Vector<Integer>> getClusters(Vector<Cluster> vtCluster){
		Vector<Vector<Integer>> Clusters = new Vector();		
		Vector<Integer> TempCluster;
		for(int count = 0; count < SimilarlityMatrixSize; count++){
			TempCluster = new Vector();
			for(int count1 = 0; count1 < SimilarlityMatrixSize ;count1++){
				if(count==NeuClusterMark[count1]) TempCluster.add(count1);
			}
			if(TempCluster.size()!=0) Clusters.add(TempCluster);
		}
		return Clusters;
	}*/
	/**
	 * The function not finished.
	 * @param para
	 * @param vtCluster
	 * @return
	 */
	public double DBI(Vector<Vector<Integer>> vtCluster){
		//Vector<InputData>
		//InputData clustercentroid = new InputData();
		Vector<InputData> ClusterCentroids = new Vector();
		listNeuron = para.getListNeuron();
		HashMap<Integer,HashMap<String,Integer>> AllCala = new HashMap();
		//HashMap<String,Integer> cala= new HashMap();
		
		Vector<Vector<Integer>> Clusters = new Vector();
		ArrayList<InputData> ListTrnInput = para.getListTrnInput();	
		for(Vector<Integer> cluster : vtCluster){
			Vector<Integer> temp = new Vector();
			for(Integer neu : cluster){
				for(int x : listNeuron.get(neu).iPrjInput){
					temp.add(x);
				}				
			}
			Clusters.add(temp); //Clusters
		}
	
		String [] Centroid = new String[para.getSAttrTypeAry().length];
		for(int j=0;j<para.getSAttrTypeAry().length;j++){
			Centroid[j] = "0";
		}
		
		for(Vector<Integer> cluster : Clusters){
			for(Integer x : cluster){
				for(int i=0;i<listNeuron.get(x).getSAnchor().length;i++){
					if (para.getSAttrTypeAry()[i] == VALUE_TYPE.integer || para.getSAttrTypeAry()[i] == VALUE_TYPE.real){
						float temp = Float.parseFloat(Centroid[i]);
						Centroid[i] = String.valueOf(temp+Float.parseFloat(ListTrnInput.get(x).getSValue()[i]));					
					}else{
						HashMap<String,Integer> cala;
						if(AllCala.get(i)!=null) cala = AllCala.get(i);
						else{
							cala= new HashMap();
							AllCala.put(i, cala);
						}
						
						String key = ListTrnInput.get(x).getSValue()[i];
						if(cala.get(key)!=null){
							int count = cala.get(key);
							cala.remove(key);
							cala.put(key, count+1);
						}else{
							cala.put(key, 1);
						}
						
						AllCala.remove(i);
						AllCala.put(i, cala);
						String maxs="";
						int maxi=0;
						for(String k : cala.keySet()){
							if(cala.get(k)>maxi){
								maxi=cala.get(k);
								maxs=k;
							}
						}
						Centroid[i] = maxs;
					}
				}
			}
			ClusterCentroids.add(new InputData(0,"",Centroid,""));			
		}
		
		/*
		for(){
			
		}
		*/			
		return 0;
	}
	/**
	 * ClusterCentroidDistanceAll
	 * @param cluster
	 * @return
	 */
	public double CCDAll(Vector<Integer> cluster){
		double distance =0;		
		return distance;
	}
	/**
	 * 
	 * @param Clusters
	 * @param Edu
	 * @return
	 */
	public double ConnIndexAllEdu(Vector<Vector<Integer>> Clusters,double[][] Edu){
		return ConnIndexHybrid(Clusters,Edu,0,0);
	}
	/**
	 * 
	 * @param Clusters
	 * @param Edu
	 * @return
	 */
	public double ConnIndex(Vector<Vector<Integer>> Clusters,double[][] Edu){
		return ConnIndexHybrid(Clusters,Edu,1,0);
	}
	/**
	 * 
	 * @param para
	 * @param Clusters
	 * @param Edu
	 * @return
	 */
	public double ConnIndexHybrid(Vector<Vector<Integer>> Clusters,double[][] Edu,double peta,double gamma){
		listNeuron = para.getListNeuron();
				
		//【【 計算指標 】】
		double Conn_Index=0;
		double Intra_Conn=0;				
		double Inter_Conn=0;
		//double peta= 0.7;
		//double gamma =0.3;

		for(Vector<Integer> vc1 : Clusters){		
			
			//Intra_Conn
			double IaCUpSec=0;
			double IaCDn=0;
			double IaEdu=0;			
			double IaCUpThi=0;
			double TriIaCDn=0;
			int Neighbors=0;		
			if(vc1.size() > 1){
				for(Integer na1:vc1){
						for(Integer na2:vc1){
							if(listNeuron.get(na1).iNeuronID!=listNeuron.get(na2).iNeuronID){
								//Second BMU
								if(listNeuron.get(na1).RFij.get(listNeuron.get(na2).iNeuronID)!=null)
									IaCUpSec+=listNeuron.get(na1).RFij.get(listNeuron.get(na2).iNeuronID);								
								//Third BMU
								if(listNeuron.get(na1).RFik.get(listNeuron.get(na2).iNeuronID)!=null){
									IaCUpThi+=listNeuron.get(na1).RFik.get(listNeuron.get(na2).iNeuronID);
								}							
								//Euclidean distance
								IaEdu+=Edu[na1][na2];
								Neighbors++;
							}
						}
						IaCDn+=listNeuron.get(na1).RFi;							
				}
			}else{			
				IaCDn=IaCUpSec=listNeuron.get(vc1.get(0)).RFi; 		//Second BMU				
				IaCUpThi=listNeuron.get(vc1.get(0)).RFi;			//Third BMU				
				IaEdu=0; Neighbors=1;								//Edu
			}
			
			if(IaCDn!=0){
				Intra_Conn += (peta*(IaCUpSec/IaCDn) + gamma*(IaCUpThi/IaCDn) + (1-peta-gamma)*(IaEdu/Neighbors));
				//Intra_Conn += IaCUpSec/IaCDn;
			}
																	
			//Inter_Conn				
			double IrC_max=0;
			double IrEdu_max=0;
			double IrC_maxk=0;
			for(Vector<Integer> vc2 : Clusters){						
				double IrCUpSec=0;
				double IrCDnSec=0;	
				
				double IrEdu=0;
				double IrCUpThi=0;
				double IrCDnThi=0;
				int IrCount =0;
				if(vc1!=vc2){ ///須確定是否有問題
					for(Integer nr1:vc1){
						double EduMin=0;					
						for(Integer nr2:vc2){							
							//判斷是否為邊界上的神經元
							boolean isnbr=false;
							int nbID[]=listNeuron.get(nr1).getINeigborID();
							for(int nb=0;nb<4;nb++){
								if(nbID[nb]==nr2) isnbr=true;
								if((nb==1 || nb==2) && nbID[nb]!=-1){
									int[] nbnbID = listNeuron.get(nbID[nb]).getINeigborID();
									if(nbnbID[3]==nr2 || nbnbID[4]==nr2) isnbr=true;
								}
							}
							
							if(isnbr){
								//Second BMU									
								if(listNeuron.get(nr1).RFij.get(listNeuron.get(nr2).iNeuronID)!=null){
									IrCUpSec+=listNeuron.get(nr1).RFij.get(listNeuron.get(nr2).iNeuronID);
									IrCDnSec+=listNeuron.get(nr1).RFij.get(listNeuron.get(nr2).iNeuronID);									
								}
								if(listNeuron.get(nr2).RFij.get(listNeuron.get(nr1).iNeuronID)!=null)
									IrCDnSec+=listNeuron.get(nr2).RFij.get(listNeuron.get(nr1).iNeuronID);
								//Third BMU
								if(listNeuron.get(nr1).RFik.get(listNeuron.get(nr2).iNeuronID)!=null){
									IrCUpThi+=listNeuron.get(nr1).RFik.get(listNeuron.get(nr2).iNeuronID);
									IrCDnThi+=listNeuron.get(nr1).RFik.get(listNeuron.get(nr2).iNeuronID);									
								}
								if(listNeuron.get(nr2).RFik.get(listNeuron.get(nr1).iNeuronID)!=null)
									IrCDnThi+=listNeuron.get(nr2).RFik.get(listNeuron.get(nr1).iNeuronID);
								//Euclidean distance
								if(EduMin < Edu[nr1][nr2]) EduMin = Edu[nr1][nr2];									
							}									
						}
						IrEdu += EduMin;
						IrCount++;
					}
				}							
				if(IrCDnSec!=0) if(IrC_max<(IrCUpSec/IrCDnSec)) IrC_max = IrCUpSec/IrCDnSec;	//Second BMU	
				if(IrCDnThi!=0) if(IrC_maxk<(IrCUpThi/IrCDnThi)) IrC_maxk = IrCUpThi/IrCDnThi;	//Third BMU
				if(IrCount!=0) IrEdu_max = IrEdu/IrCount;										//Euclidean distance
			}									
			Inter_Conn += (peta*IrC_max + gamma*IrC_maxk + (1-peta-gamma)*IrEdu_max);
			//Inter_Conn += IrC_max;
		}				
		 
		System.out.println("Ia:"+Intra_Conn+" Ir:"+Inter_Conn);
		Intra_Conn = Intra_Conn / (float)Clusters.size();
		Inter_Conn = Inter_Conn / (float)Clusters.size();
		Conn_Index = Intra_Conn * (1 - Inter_Conn);									
		return Conn_Index;
	}
	
	/**
	 * 
	 * @param para
	 * @param vtCluster
	 * @return
	 */
	public double Silhouette(Vector<Vector<Integer>> vtCluster){
		double Sil = 0;
		Comm comm = new Comm();
		listNeuron = para.getListNeuron(); //listNeuron
		int totalsample = para.getITotSample();
		
		
		//Vector<Vector<Integer>> Clusters = vtCluster;
		
		Vector<Vector<Integer>> Clusters = new Vector<Vector<Integer>>();
		ArrayList<InputData> ListTrnInput = para.getListTrnInput();	
		
		for(Vector<Integer> cluster : vtCluster){ 			//取各群
			Vector<Integer> temp = new Vector<Integer>();
			System.out.println("Cluster: " + cluster);
			for(Integer neu : cluster){ 					//取得群內神經元編號
				for(int x : listNeuron.get(neu).iPrjInput){ //取得神經元內的資料編號
					temp.add(x);
				}				
			}
			Clusters.add(temp); //Clusters
		}
		
		
		
		for(Vector<Integer> vt : Clusters){							
			for(int xi : vt){
				double b = 0;
				double a = 0;
				
				//到群內
				for(int xj : vt){
					if(xi != xj){
						a += comm.getEuclideanDistance(ListTrnInput.get(xi), ListTrnInput.get(xj), para) / ListTrnInput.get(xj).getSValue().length;
						//a += comm.getEuclideanDistance(listNeuron.get(xi), listNeuron.get(xj), para);
					}
				}
				//System.out.println("a/size: " + a +"/" + vt.size() + "   " + (a/vt.size()));
				a = a/vt.size();
								
				//到其他群
				double min = Double.POSITIVE_INFINITY;				
				for(Vector<Integer> TempCluster : Clusters){
					double TempValue = 0.0;
					if(vt != TempCluster){
						for(int xj : TempCluster){
							//System.out.println(comm.getEuclideanDistance(ListTrnInput.get(xi), ListTrnInput.get(xj), para));
							TempValue += comm.getEuclideanDistance(ListTrnInput.get(xi), ListTrnInput.get(xj), para)/ListTrnInput.get(xj).getSValue().length;
							//cl += comm.getEuclideanDistance(listNeuron.get(xi), listNeuron.get(xj), para);
						}
						TempValue = TempValue/TempCluster.size();
						if(min > TempValue) min = TempValue;
					}
				}				
				b = min;
				
				double max =0.0;
				if(a>b) max = a;
				else max = b;
				//System.out.println("a:" + a + "  b:" + b +"  Sil : " + ((b-a)/max));
				//System.out.println("Sil : " + ((b-a)/max));
				if (max < Double.POSITIVE_INFINITY) Sil += (b-a)/max;
			}			
		}
		
		//System.out.println("Sil:"+Sil+"  Cluste size:" + Clusters.size());
		return Sil/totalsample;
	}
}
