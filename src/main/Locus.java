package src.main;

public enum Locus {
	A{
		public int[][] getPockets(){
			return new int[][] {
					{63, 66, 99, 163, 167, 171}, // pocket A
					{9, 63, 66, 67, 70, 99},     // pocket B
					{9, 70, 74, 97},             // pocket C
					{99, 114, 155, 156},         // pocket D
					{97, 114, 152, 156},         // pocket E
					{77, 80, 116}                // pocket F
			} ;
		}
		public String getRefSeq(){
			return "01:01:01:01";
		}
	}, 
	B{
		public int[][] getPockets(){
			return new int[][] {
					{63, 66, 99, 163, 167, 171}, // pocket A
					{9, 63, 66, 67, 70, 99},     // pocket B
					{9, 70, 74, 97},             // pocket C
					{99, 114, 155, 156},         // pocket D
					{97, 114, 152, 156},         // pocket E
					{77, 80, 116}                // pocket F
			};
		}
		public String getRefSeq(){
			return "07:02:01";
		}
	},  
	Cw{
		public int[][] getPockets(){
			return new int[][] {
					{63, 66, 99, 163, 167, 171}, // pocket A
					{9, 63, 66, 67, 70, 99},     // pocket B
					{9, 70, 74, 97},             // pocket C
					{99, 114, 155, 156},         // pocket D
					{97, 114, 152, 156},         // pocket E
					{77, 80, 116}                // pocket F
			};
		}
		public String getRefSeq(){
			return "01:02:01";
		}
	},
    C{
        public int[][] getPockets(){
            return new int[][] {
                    {63, 66, 99, 163, 167, 171}, // pocket A
                    {9, 63, 66, 67, 70, 99},     // pocket B
                    {9, 70, 74, 97},             // pocket C
                    {99, 114, 155, 156},         // pocket D
                    {97, 114, 152, 156},         // pocket E
                    {77, 80, 116}                // pocket F
            };
        }
        public String getRefSeq(){
            return "01:02:01";
        }
    },
	DRB1{
		public int[][] getPockets(){
			return new int[][] {
					{82,85,86,89},          // pocket 1
					{},                  //
					{},                  //
					{13,26,70,71,74,78}, // pocket 4
					{},                  //
					{11},           // pocket 6
					{28,47,61,67,71},    // pocket 7
					{},                  //
					{9,57,60,61}      // pocket 9
			};
		}
		public String getRefSeq(){
			return "01:01:01";
		}
	},  
	DQB1{
		public int[][] getPockets(){
			return new int[][] {
                    {82,85,86,89},          // pocket 1
                    {},                  //
                    {},                  //
                    {13,26,70,71,74,78}, // pocket 4
                    {},                  //
                    {11},           // pocket 6
                    {28,47,61,67,71},    // pocket 7
                    {},                  //
                    {9,57,60,61}      // pocket 9
			};
		}
		public String getRefSeq(){
			return "05:01:01:01";
		}
	},  
	DPB1{
		public int[][] getPockets(){
			return new int[][] {
					{87,84},             // pocket 1
					{},                  //
					{},                  //
					{13,69,76,68,72,24}, // pocket 4
					{},                  //
					{9,11,28},           // pocket 6
					{26,59,69,45,65},    // pocket 7
					{},                  //
					{9,58,55,35,36}      // pocket 9
			};
		}
		public String getRefSeq(){
			return "01:01:01";
		}
	},
	DQA1{
		public int[][] getPockets(){
			return new int[][] {
					{34,46,56,35,55,57,27},  // pocket 1
					{},                      //
					{},                      //
					{11,65,14},              // pocket 4
					{},                      //
					{14,68,66,69},           // pocket 6
					{68,72},                 // pocket 7
					{},                      //
					{75,76,72,79}            // pocket 9
			};
		}
		public String getRefSeq(){
			return "01:01:01";
		}
	},
	DPA1{
		public int[][] getPockets(){
			return new int[][] {
					{31,43,53,32,52,54,24},  // pocket 1
					{},                      //
					{},                      //
					{9,62,11},               // pocket 4
					{},                      //
					{11,65,62,66},           // pocket 6
					{65,69},                 // pocket 7
					{},                      //
					{72,73,69,76}            // pocket 9
			};
		}
		public String getRefSeq(){
			return "01:03:01";
		}
	},
	MICA{
		public int[][] getPockets(){
			return new int[][] {{}};
		}
		public String getRefSeq(){
			return "001";
		}
	},  
	MICB{
		public int[][] getPockets(){
			return new int[][] {{}};
		}
		public String getRefSeq(){
			return "001";
		}
	},  
	ERR{
		public int[][] getPockets(){
			return new int[][] {{}};
		}
		public String getRefSeq(){
			return "";
		}
	},
    DMA1{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    DMB1{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    DOA{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01";
        }
    },
    DOB{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    DRA{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    E{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    F{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    G{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    TAP1{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01:01";
        }
    },
    TAP2{
        public int[][] getPockets(){
            return new int[][] {{}};
        }
        public String getRefSeq(){
            return "01:01:01";
        }
    };
	
	public abstract int[][] getPockets();
	public abstract String getRefSeq();
	
};
