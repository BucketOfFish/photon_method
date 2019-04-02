
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip> 

vector< vector<int> > GetNThingsFromAGroup(int n, int m) {
	// we will take n objects from a group of m objects
	vector< vector<int> > list;
	if (n==1) {
		for (int i0=0;i0<m;i0++) {
			list.push_back(vector <int> ());
			int current_list = list.size()-1;
			list.at(current_list).push_back(i0);
		}
	}
	else if (n==2) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				list.push_back(vector <int> ());
				int current_list = list.size()-1;
				list.at(current_list).push_back(i0);
				list.at(current_list).push_back(i1);
			}
		}
	}
	else if (n==3) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					list.push_back(vector <int> ());
					int current_list = list.size()-1;
					list.at(current_list).push_back(i0);
					list.at(current_list).push_back(i1);
					list.at(current_list).push_back(i2);
				}
			}
		}
	}
	else if (n==4) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						list.push_back(vector <int> ());
						int current_list = list.size()-1;
						list.at(current_list).push_back(i0);
						list.at(current_list).push_back(i1);
						list.at(current_list).push_back(i2);
						list.at(current_list).push_back(i3);
					}
				}
			}
		}
	}
	else if (n==5) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						for (int i4=i3+1;i4<m;i4++) {
							list.push_back(vector <int> ());
							int current_list = list.size()-1;
							list.at(current_list).push_back(i0);
							list.at(current_list).push_back(i1);
							list.at(current_list).push_back(i2);
							list.at(current_list).push_back(i3);
							list.at(current_list).push_back(i4);
						}
					}
				}
			}
		}
	}
	else if (n==6) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						for (int i4=i3+1;i4<m;i4++) {
							for (int i5=i4+1;i5<m;i5++) {
								list.push_back(vector <int> ());
								int current_list = list.size()-1;
								list.at(current_list).push_back(i0);
								list.at(current_list).push_back(i1);
								list.at(current_list).push_back(i2);
								list.at(current_list).push_back(i3);
								list.at(current_list).push_back(i4);
								list.at(current_list).push_back(i5);
							}
						}
					}
				}
			}
		}
	}
	else if (n==7) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						for (int i4=i3+1;i4<m;i4++) {
							for (int i5=i4+1;i5<m;i5++) {
								for (int i6=i5+1;i6<m;i6++) {
									list.push_back(vector <int> ());
									int current_list = list.size()-1;
									list.at(current_list).push_back(i0);
									list.at(current_list).push_back(i1);
									list.at(current_list).push_back(i2);
									list.at(current_list).push_back(i3);
									list.at(current_list).push_back(i4);
									list.at(current_list).push_back(i5);
									list.at(current_list).push_back(i6);
								}
							}
						}
					}
				}
			}
		}
	}
	else if (n==8) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						for (int i4=i3+1;i4<m;i4++) {
							for (int i5=i4+1;i5<m;i5++) {
								for (int i6=i5+1;i6<m;i6++) {
									for (int i7=i6+1;i7<m;i7++) {
										list.push_back(vector <int> ());
										int current_list = list.size()-1;
										list.at(current_list).push_back(i0);
										list.at(current_list).push_back(i1);
										list.at(current_list).push_back(i2);
										list.at(current_list).push_back(i3);
										list.at(current_list).push_back(i4);
										list.at(current_list).push_back(i5);
										list.at(current_list).push_back(i6);
										list.at(current_list).push_back(i7);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else if (n==9) {
		for (int i0=0;i0<m;i0++) {
			for (int i1=i0+1;i1<m;i1++) {
				for (int i2=i1+1;i2<m;i2++) {
					for (int i3=i2+1;i3<m;i3++) {
						for (int i4=i3+1;i4<m;i4++) {
							for (int i5=i4+1;i5<m;i5++) {
								for (int i6=i5+1;i6<m;i6++) {
									for (int i7=i6+1;i7<m;i7++) {
										for (int i8=i7+1;i8<m;i8++) {
											list.push_back(vector <int> ());
											int current_list = list.size()-1;
											list.at(current_list).push_back(i0);
											list.at(current_list).push_back(i1);
											list.at(current_list).push_back(i2);
											list.at(current_list).push_back(i3);
											list.at(current_list).push_back(i4);
											list.at(current_list).push_back(i5);
											list.at(current_list).push_back(i6);
											list.at(current_list).push_back(i7);
											list.at(current_list).push_back(i8);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//for (int i=0;i<list.size();i++) {
	//	for (int j=0;j<list.at(i).size();j++) {
	//		std::cout << list.at(i).at(j) << ", ";
	//	}
	//	std::cout << std::endl;
	//} 
	return list;
}
