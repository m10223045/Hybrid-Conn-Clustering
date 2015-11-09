/**
 * 
 */
package src;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JOptionPane;

import src.Constant.UI_SET;

public class ClusterCONN {
	//private FileSystem filesystemA;
	//private FileSystem filesystemB;
	//private FileSystem filesystemC;
	private Vector<Cluster> vtTempCluster;
	private Vector<Double> vtConn_Index;
	private Vector<Cluster> vtCluster;
	private int iClusterNeed;
	private ArrayList<Neuron> listNeuron = null;
	private ExprParm para = null;
	/**
	 * ClusterCONN Constructor
	 * @param para
	 * @param vtCluster
	 * @param iClusterNeed
	 */
	public ClusterCONN(ExprParm para, Vector<Cluster> vtCluster, int iClusterNeed) {
		try {
			//filesystemA = new FileSystem("D:\\", "SimilarityMartiex", "csv");
			//filesystemB = new FileSystem("D:\\", "CONNcluster", "csv");
			//filesystemC = new FileSystem("D:\\", "Check", "csv");

			vtTempCluster = new Vector<Cluster>();
			vtConn_Index = new Vector<Double>();

			this.listNeuron = para.getListNeuron();
			this.para = para;
			this.vtCluster = vtCluster;
			this.iClusterNeed = iClusterNeed;
		} catch (Exception e) {
		}
	}
	
	/**
	 * Merger two cluster to one
	 * @param cluster1
	 * @param cluster2
	 */
	private void MergedCluster(Cluster cluster1, Cluster cluster2) {
		Cluster master = (cluster1.iClusterID < cluster2.iClusterID) ? cluster1
				: cluster2;
		Cluster slave = (cluster1.iClusterID < cluster2.iClusterID) ? cluster2
				: cluster1;
		Iterator<Neuron> itrNeuron = slave.vtNeuron.iterator();
		while (itrNeuron.hasNext())
			master.addNeuron(itrNeuron.next());
		slave.clearNeuron(); // clear neuron ID
	}
	/**
	 * Here using Array to mark the cluster id for each neurons.  In the finally, we choose the best cluster form vtConn_index then building the Cluster object. 
	 * @throws IOException
	 */
	public void doClustering() throws IOException {
		//String[] Control = JOptionPane.showInputDialog("�п�J����Ѽƣ\/�] ex:0.7/0.2").split("/");
		
		ClusterValidityIndices CVI = new ClusterValidityIndices(para); //���s���T�ʫ��Фu��
		Vector<Cluster> vtTempCluster = new Vector<Cluster>(); // �s��Ȯɸs��
		Vector<Double> vtConn_Index = new Vector<Double>();
		//Vector<Double> vtSilhouette = new Vector<Double>();
		//Vector<Double> vtConnIndexHyban = new Vector<Double>();
		int SimilarlityMatrixSize = listNeuron.size(); 
		int TrueClusterSize;
		int CurrentClusterNumber = SimilarlityMatrixSize;

		// �p��RFij,RFi
		for (Neuron Neuron_i : listNeuron) {
			System.out.println("Neurn " + Neuron_i.iNeuronID);
			System.out.println("Input:" + Neuron_i.iPrjInput.size()
					+ " Second: " + Neuron_i.iPrjInputSecondBMUID.size()
					+ " Third:" + Neuron_i.iPrjInputTirdBMUID.size());
			
			Neuron_i.RFi = Neuron_i.iPrjInputSecondBMUID.size(); // set RFi that equal to the number of input data	
			
			// �p��RFij
			for (Neuron Neuron_j : listNeuron) {
				int sum = 0;
				for (Integer i : Neuron_i.iPrjInputSecondBMUID) {
					if (i == Neuron_j.iNeuronID)
						sum += 1;
				}
				if (sum > 0)
					Neuron_i.RFij.put(Neuron_j.iNeuronID, sum);
			}
			System.out.println("");
			
			// �p��RFik
			for (Neuron Neuron_k : listNeuron) { 
				int sumK = 0;
				for (Integer i : Neuron_i.iPrjInputTirdBMUID) {
					if (i == Neuron_k.iNeuronID)
						sumK += 1;
				}
				if (sumK > 0)
					Neuron_i.RFik.put(Neuron_k.iNeuronID, sumK);
			}
		}
		
		System.out.println("OriginalClusterSize = " + SimilarlityMatrixSize);
		
		// �ŧi�ۦ��ʯx�}
		double[][] S = new double[SimilarlityMatrixSize][SimilarlityMatrixSize]; 
		double[][] EuclideanDistance = new double[SimilarlityMatrixSize][SimilarlityMatrixSize];
		double ConnMax = 0;
		
		//�ŧi����Ѽ� �\(alpha) and �](beta), �_(delta) and �^(gamma) �Ψӱ��� Cluster Validity Indices   
		//double alpha = Double.parseDouble(Control[0]); 	// Second 
		//double beta = Double.parseDouble(Control[1]); 	// Third 
		double alpha = Double.parseDouble(para.getFrmParent().getUIsetStr(UI_SET.Alpha)); 	// Second BMU Control Parameter
		double beta = Double.parseDouble(para.getFrmParent().getUIsetStr(UI_SET.Beta)); 	// Third BMU Control Parameter
		double delta = alpha;
		double gamma = beta;
		
		if((alpha+beta)>1){
			JOptionPane.showMessageDialog(null, "Error: ����Ѽ� �\+�] ���i  >1", "ConnHybrid Error", JOptionPane.ERROR_MESSAGE );
		}else{
			for (int a = 0; a < SimilarlityMatrixSize; a++) {
				for (int b = a; b < SimilarlityMatrixSize; b++) {
								
					//Euclidean Distance Matrix
					if (a == b) EuclideanDistance[a][b] = Double.POSITIVE_INFINITY;
					else EuclideanDistance[a][b] = EuclideanDistance[b][a] = 1-Comm.getEuclideanDistance(listNeuron.get(a), listNeuron.get(b), para)/listNeuron.get(a).getSAnchor().length;
					
					
					int rfij = 0;
					int rfji = 0;
					int rfik = 0;
					int rfki = 0;
					if (a != b /* && listNeuron.get(a).RFi != 0 */) {
											
						//Second BMU
						if (listNeuron.get(a).RFij.get(listNeuron.get(b).iNeuronID) != null)
							rfij = listNeuron.get(a).RFij.get(listNeuron.get(b).iNeuronID);
						if (listNeuron.get(b).RFij.get(listNeuron.get(a).iNeuronID) != null)
							rfji = listNeuron.get(b).RFij.get(listNeuron.get(a).iNeuronID);
						double connNorUp = (double) (rfij + rfji);
						double connNorDown = (double) para.getITotSample();
						double connNor = 0.0;
						if (ConnMax < connNorUp) ConnMax = connNorUp;
						if (connNorDown != 0) connNor = connNorUp / connNorDown;


						// Third BMU
						if (listNeuron.get(a).RFik.get(listNeuron.get(b).iNeuronID) != null)
							rfik = listNeuron.get(a).RFik.get(listNeuron.get(b).iNeuronID);
						if (listNeuron.get(b).RFik.get(listNeuron.get(a).iNeuronID) != null)
							rfki = listNeuron.get(b).RFik.get(listNeuron.get(a).iNeuronID);
						double connNorUpTri = (double) (rfik + rfki);
						double connNorTri = 0.0;
						if (connNorDown != 0) {
							connNorTri = connNorUpTri / connNorDown;
							connNor = connNorUp / connNorDown;
						}
						
						/* Original Method */
						// S[a][b] = S[b][a] = rfij + rfji; 
						
						S[a][b] = S[b][a] = alpha * connNor + beta * connNorTri + (1 - beta - alpha) * EuclideanDistance[a][b];
					} else {
						S[a][b] = S[b][a] = -1;
					}
				}
				System.out.println("");
				//System.out.println("Length:"+ listNeuron.get(a).getSAnchor().length);
			}

			
			int[] ClusterMark = new int[SimilarlityMatrixSize];				// �ϥΤ@���}�C�����U��Neuron�����@�s
			Vector<int[]> ClusterMarkForEachPhase = new Vector<int[]>();	// �ΨӰO���C���q���s���p

			// Neuron���ݸs��l�Ƴ]�w�A�]�w�U�Ӧۤv���@�s�A�h����Neuron(�аO��-1)
			for (int nc = 0; nc < SimilarlityMatrixSize; nc++) {
				if (listNeuron.get(nc).iPrjInput.size() != 0) {
					ClusterMark[nc] = nc;
				} else {
					CurrentClusterNumber--;
					ClusterMark[nc] = -1;
				}
			}
			System.out.println("Cleard Neurons:" + CurrentClusterNumber);
			TrueClusterSize = CurrentClusterNumber;

			/* START CLUSTERING */
			while (CurrentClusterNumber > 1) {
				boolean FindMaxtage = false;
				int i = 0, j = 0;
				double max = -1;

				System.out.println("CurrentClusterNumber:" + CurrentClusterNumber);
				//filesystemA.writeSimilraityMatrix(S, SimilarlityMatrixSize, CurrentClusterNumber);

				// Search Max Value
				for (int sRow = 0; sRow < SimilarlityMatrixSize; sRow++) {
					if (ClusterMark[sRow] != -1) {
						for (int sCol = sRow + 1; sCol < SimilarlityMatrixSize; sCol++) {
							if (S[sRow][sCol] > max && sRow != sCol
									&& ClusterMark[sCol] > 0) { // != -1 > 0
								i = sRow;
								j = sCol;
								max = S[sRow][sCol];
								FindMaxtage = true;
							}
						}
					}
				}

				//filesystemB.writeLine("MAX:" + max + "         " + i + " Mager " + j);			
				if (FindMaxtage) {
					// �p���s�ƶq
					double FirstClusterSize = 0, SecondClusterSize = 0;
					int count = 0;

					while (count < SimilarlityMatrixSize) {
						if (ClusterMark[i] == ClusterMark[count])
							FirstClusterSize++;
						if (ClusterMark[j] == ClusterMark[count])
							SecondClusterSize++;
						count++;
					}

					// ��s�X�֫�Neu�����s�����A�N�s���j���M�֤J�s���p��
					count = 0;
					while (count < SimilarlityMatrixSize) {
						if (i > j) {
							if (ClusterMark[count] == i)
								ClusterMark[count] = j;
						} else {
							if (ClusterMark[count] == j)
								ClusterMark[count] = i;
						}
						// System.out.print(ClusterMark[count]+"  ");
						count++;
					}
					ClusterMarkForEachPhase.add(Arrays.copyOf(ClusterMark, ClusterMark.length)); // ������e���s��T

					// ��s�x�}
					int k = 0;
					if (i > j) k = j;
					else k = i;
					for (int t = 0; t < SimilarlityMatrixSize; t++) {
						if ((t != i && t != j) && t != k && ClusterMark[t] != -1) {
							S[k][t] = S[t][k] = ((FirstClusterSize / (FirstClusterSize + SecondClusterSize))
									* S[i][t] + (SecondClusterSize / (FirstClusterSize + SecondClusterSize))
									* S[j][t]);
						}
					}
					for (int d = 0; d < SimilarlityMatrixSize; d++) {
						if (k == i) S[d][j] = S[j][d] = -1;
						else if (k == j) S[d][i] = S[i][d] = -1;
					}

					// �Ȯɹ갵��e���s
					Vector<Vector<Integer>> Clusters = new Vector();
					Vector<Integer> TempCluster;
					for (count = 0; count < SimilarlityMatrixSize; count++) {
						TempCluster = new Vector<Integer>();
						for (int count1 = 0; count1 < SimilarlityMatrixSize; count1++) {
							if (count == ClusterMark[count1])
								TempCluster.add(count1);
						}
						if (TempCluster.size() != 0)
							Clusters.add(TempCluster);
					}
					CurrentClusterNumber = Clusters.size();
					// filesystemB.writeLine("Current Number of Clusters=" + CurrentClusterNumber + " ArrayIndex=" + (ClusterMarkForEachPhase.size() - 1) + "  AllClusterSize=" + Clusters.size());

					// �i�i �p����� �j�j
					vtConn_Index.add(CVI.ConnIndexHybrid(Clusters, EuclideanDistance, delta, gamma));
					// vtConn_Index.add(CVI.ConnIndexHyban(Clusters, Edu, delta, gamma));
					// vtSilhouette.add(CVI.Silhouette(Clusters));
					// vtConnIndexHyban.add(CVI.ConnIndexHyban(Clusters, Edu, delta, gamma));
				}
			}

			System.out.println("");
			System.out.println("Total Neuron=" + SimilarlityMatrixSize
					+ " **** NeedClusterNumber=" + iClusterNeed
					+ " TrueClusterSize=" + TrueClusterSize);

			// �P�_�O�_�ϥΪ̬O�_���]�w���s�ƶq
			int better = 0;
			if (iClusterNeed != 0 && iClusterNeed <= TrueClusterSize) better = vtConn_Index.size() - iClusterNeed;
			else{
				// Find Array Index of Better Clustering 
				for (int bt = 0; bt < vtConn_Index.size() - 1; bt++) {
					if (vtConn_Index.get(better) < vtConn_Index.get(bt)) better = bt;
					
					//filesystemB.writeLine((vtConn_Index.size() - bt) + "," + vtConn_Index.get(bt) /*+ "," + vtSilhouette.get(bt)+ "," + vtConnIndexHyban.get(bt)*/);
					System.out.print(vtConn_Index.get(bt) + ","); // �L�X�Ҧ�����
				}			
			}
			System.out.println("vtC_Neu:" + ClusterMarkForEachPhase.size() + " **** Conn_list Size=" + vtConn_Index.size() + " and Better in " + better);

			// ��l�ƼȮɸs��
			int number = 0;
			for (Neuron neu : listNeuron) {
				vtTempCluster.add(new Cluster(neu, number++, para));
			}

			// �̩ҳ̨Τ��s�i��s���X�֡A-1��ܤw�X�֧���
			int[] finalNeuCluster = ClusterMarkForEachPhase.get(better);
			for (int c = 0; c < SimilarlityMatrixSize; c++) {
				if (finalNeuCluster[c] != -1) {
					for (int search = c + 1; search < SimilarlityMatrixSize; search++) {
						if (finalNeuCluster[search] != -1
								&& finalNeuCluster[c] == finalNeuCluster[search]) {
							MergedCluster(vtTempCluster.get(c),
									vtTempCluster.get(search));
							finalNeuCluster[search] = -1;
						}
					}
				}
				finalNeuCluster[c] = -1;
			}

			// ��ڸs���]�w�A�аO��Neuron��-1
			for (Cluster Ct : vtTempCluster) {
				if (Ct.vtNeuron.size() > 1) vtCluster.add(Ct);
				else if (Ct.vtNeuron.size() == 1) {
					if (Ct.vtNeuron.get(0).iPrjInput.size() > 0) vtCluster.add(Ct);
					else Ct.vtNeuron.get(0).iClusterID = -1;
				}
			}

			for (int qq = 0; qq < para.getSClassAry().length; qq++)
				System.out.println("index:" + qq + "  " + para.getSClassAry()[qq]);

			System.out.println("Entropy:" + para.getPerformance().getClusteringEntropy(para, vtCluster));		
		}//if end 

		//System.out.println("Silhouette:" +
		//para.getPerformance().getSilhouette(para, vtCluster));
		//return 0;
	}

}
