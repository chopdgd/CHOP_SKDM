package src.gui;

import src.main.Run;
import src.utils.Settings;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;
import java.util.Calendar;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;


public class Window extends JPanel implements ActionListener {
	private static final long serialVersionUID = -2834115717777574977L;
    protected static String[] tfTxts;
    protected static JTextField[] textFields;
    protected static final String btRun = "buttonRun";
    protected static JTextArea hiddenPane;
    protected static JTextArea hiddenDetailsPane;
    protected static JTextArea inputPane1;
    protected static JTextArea inputPane2;
    protected JLabel actionLabel;
    protected static JTabbedPane tabbedPane;

    public Window() {
        setLayout(new BorderLayout());

        // RUN button
        JButton button = new JButton();
        button.setText("     R U N  . . .");
        button.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        button.setToolTipText("Run current analysis on a new tab, or last analysis if input is blank");
        button.setMargin(new Insets(0,250,0,250));
        button.setActionCommand(btRun);
        button.addActionListener(this);

        // labels
        JLabel[] textFieldLabels = new JLabel[tfTxts.length];
        for (int i = 0; i < textFieldLabels.length; i++){
        	textFields[i].setActionCommand(tfTxts[i]);
        	textFields[i].addActionListener(this);
        	textFieldLabels[i] = new JLabel(tfTxts[i] + ": ");
        	textFieldLabels[i].setLabelFor(textFields[i]);
        }

        // action event label
        actionLabel = new JLabel("");
        actionLabel.setBorder(BorderFactory.createEmptyBorder(5,0,5,0));

        // lay out text controls and labels
        JPanel textControlsPane = new JPanel();
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        textControlsPane.setLayout(gridbag);

        addLabelTextRows(textFieldLabels, textFields, gridbag, textControlsPane);

        c.gridwidth = GridBagConstraints.REMAINDER; //last
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 1.0;
        textControlsPane.add(actionLabel, c);
        textControlsPane.add(button, c);
        textControlsPane.setBorder(
                BorderFactory.createCompoundBorder(
                                BorderFactory.createTitledBorder("Run"),
                                BorderFactory.createEmptyBorder(5,5,5,5)));

        // Output
        tabbedPane = new JTabbedPane();
        tabbedPane.setBorder(BorderFactory.createCompoundBorder(
        		BorderFactory.createCompoundBorder(
        				BorderFactory.createTitledBorder("Output"),
        				BorderFactory.createEmptyBorder(5,5,5,5)),
        				tabbedPane.getBorder()));
        
        // input 1
        inputPane1 = createTextArea("Copy-Paste Dataset 1 (CASES)");
        JScrollPane inputPane1Scroll = new JScrollPane(inputPane1);
        inputPane1Scroll.setPreferredSize(new Dimension(300, 400));
        inputPane1Scroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
        inputPane1Scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        
        // input 2
        inputPane2 = createTextArea("Copy-Paste Dataset 2 (CONTROLS)");
        JScrollPane inputPane2Scroll = new JScrollPane(inputPane2);
        inputPane2Scroll.setPreferredSize(new Dimension(300, 400));
        inputPane2Scroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
        inputPane2Scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        
        // hidden
        hiddenPane = createTextArea("");
        hiddenDetailsPane = createTextArea("");
        initHTML();
        
        // horizontal Line
        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, inputPane1Scroll, inputPane2Scroll);
        splitPane.setOneTouchExpandable(true);
        splitPane.setResizeWeight(0.5);
        
        // left Pane
        JPanel leftPane = new JPanel(new GridLayout(1,0));
        leftPane.add(splitPane);
        leftPane.setBorder(BorderFactory.createCompoundBorder(
                        BorderFactory.createTitledBorder("Input"),
                        BorderFactory.createEmptyBorder(5,5,5,5)));

        // right Pane
        JPanel rightPane = new JPanel(new BorderLayout());
        rightPane.add(textControlsPane, BorderLayout.PAGE_START);
        rightPane.add(tabbedPane, BorderLayout.CENTER);

        // veritcal Line
        JSplitPane splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftPane, rightPane);
        splitPane.setOneTouchExpandable(true);
        splitPane.setResizeWeight(0.5);
        
        add(splitPane2, BorderLayout.CENTER);// .LINE_START);
        //add(rightPane, BorderLayout.LINE_END);
    }

    protected JComponent makeTextPanel(String text) {
    	JEditorPane filler = new JEditorPane();
    	filler.setEditable(false);
    	filler.setContentType("text/html");
    	JScrollPane scrollPane = new JScrollPane(filler);
    	scrollPane.setPreferredSize(new Dimension(700, 900));
    	scrollPane.setMinimumSize(new Dimension(10, 10));

        JPanel panel = new JPanel(false);
        filler.setText(text);
        panel.setLayout(new GridLayout(1, 1));
        panel.add(scrollPane);
        return panel;
    }
    
    //Returns an ImageIcon, or null if the path was invalid.
    protected static ImageIcon createImageIcon(String path) {
        java.net.URL imgURL = Window.class.getResource(path);
        if (imgURL != null) {
            return new ImageIcon(imgURL);
        } else {
            System.err.println("Couldn't find file: " + path);
            return null;
        }
    }
    
    private void addLabelTextRows(JLabel[] labels,
                                  JTextField[] textFields,
                                  GridBagLayout gridbag,
                                  Container container) {
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.EAST;
        int numLabels = labels.length;

        for (int i = 0; i < numLabels; i++) {
            c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
            c.fill = GridBagConstraints.NONE;      //reset to default
            c.weightx = 0.0;                       //reset to default
            container.add(labels[i], c);

            c.gridwidth = GridBagConstraints.REMAINDER;     //end row
            c.fill = GridBagConstraints.HORIZONTAL;
            c.weightx = 1.0;
            container.add(textFields[i], c);
        }
    }

    public void actionPerformed(ActionEvent e) {
        for (int i = 0; i < tfTxts.length; i++){
        	if (tfTxts[i].equals(e.getActionCommand())){
        		JTextField source = (JTextField)e.getSource();
                actionLabel.setText(tfTxts[i] + " now set to \"" + source.getText() + "\"");
                saveSettings();
        	}
        }

        if (btRun.equals(e.getActionCommand())) {
        	saveSettings();

        	initHTML();
            String in1 = inputPane1.getText();
            String in2 = inputPane2.getText();
            if (in1.length() > 1 && in2.length() > 1)
            	Run.main(new String[] {in1, in2});
            else if (in1.length() > 1 && in2.length() < 1)
            	Run.main(new String[] {in1, null});
            else if (in1.length() < 1 && in2.length() > 1)
            	Run.main(new String[] {null, in2});
            else
            	Run.main(null);
            
            appendOutputTxt("</pre><hr><br><br></body></html>");
            appendDetailsTxt("</pre><hr><br><br></body></html>");
            
            JComponent panel1 = makeTextPanel(hiddenPane.getText());
            tabbedPane.addTab("Results", panel1);
            tabbedPane.setToolTipTextAt(tabbedPane.getTabCount()-1, "Analysis results");
            String tbtText = hiddenDetailsPane.getText();
            if (tbtText.length() > 300){
	            JComponent panel2 = makeTextPanel(tbtText);
	            tabbedPane.addTab("2x2", panel2);
	            tabbedPane.setToolTipTextAt(tabbedPane.getTabCount()-1, "Two by two contingency tables used in the analysis");
	            tabbedPane.setSelectedIndex(tabbedPane.getTabCount()-2);
            }
            else
            	tabbedPane.setSelectedIndex(tabbedPane.getTabCount()-1);
            
            actionLabel.setText("Analysis completed on " + Calendar.getInstance().getTime().toString() + ".");
        }
    }

    private final void saveSettings(){
    	LinkedList<String> ll = new LinkedList<String>();
    	for (int i = 0; i < textFields.length; i++)
    		ll.add(textFields[i].getText());
    	Settings.setSettings(ll.iterator());
    }
    
    private final void initHTML(){
    	hiddenPane.setText("<html><body style='font-family:Courier New; font-size:9px;'><pre>");
    	hiddenDetailsPane.setText("<html><body style='font-family:Courier New; font-size:9px;'><pre>");
    }
    
    private JTextArea createTextArea(String ttt){
    	JTextArea textArea = new JTextArea("");
    	textArea.setFont(new Font("Courier New", Font.PLAIN, 12));
    	textArea.setLineWrap(false);
    	textArea.setToolTipText(ttt);
    	//textArea.setPreferredSize(new Dimension(250, 250));
    	return textArea;
    }
    
    public static void appendOutputTxt(String s){
    	hiddenPane.setText(hiddenPane.getText()+s);
    }
    
    public static void appendDetailsTxt(String s){
    	hiddenDetailsPane.setText(hiddenDetailsPane.getText()+s);
    }
    
    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("SKDM HLA Tools " + Settings.VERSION);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocationByPlatform(true);
        System.out.println("JFrame: "+frame);
        try{
        	URL url = Window.class.getResource("hla_icon.GIF");
            System.out.println("url is: " + url + "^^");
        	Image image = Toolkit.getDefaultToolkit().getImage(url);
        	frame.setIconImage(image);
        }catch(Exception e){
            e.printStackTrace();
        }
        
        //settings ********************************************************************************************
		Iterator<Entry<String,Object>> it = Settings.getSettings();
		Entry<String,Object> e;
		LinkedList<String> tfTxt = new LinkedList<String>();
		LinkedList<JTextField> textField = new LinkedList<JTextField>();
		while (it.hasNext()){
			e = it.next();
			tfTxt.add(e.getKey().toString().split("!")[0].split(":")[1]);
			JTextField jtf = new JTextField(10);
			jtf.setText(e.getValue().toString());
			jtf.setToolTipText(e.getKey().toString().split("!")[1].trim());
			textField.add(jtf);
		}
		tfTxts = tfTxt.toArray(new String[tfTxt.size()]);
		textFields = textField.toArray(new JTextField[textField.size()]);
		
        //Create and set up the content pane.
        JComponent newContentPane = new Window();
        newContentPane.setOpaque(true);
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
            	Settings.load();
                createAndShowGUI();
            }
        });
    }
}
