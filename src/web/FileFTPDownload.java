package src.web;

/**
 * Created by jayaramanp on 2/4/15.
 */

import src.main.Locus;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import java.io.*;

    /**
     * A program demonstrates how to download files from remote FTP server to a local computer using Apache Commons Net API.
     * @author www.codejava.net
     */
    public class FileFTPDownload {

        public void downloadFile(Locus locus) {
            String server = "ftp.ebi.ac.uk";
            int port = 21;
            String user = "";
            String pass = "";

            FTPClient ftpClient = new FTPClient();
            try {

                ftpClient.connect(server, port);
                ftpClient.login(user, pass);
                ftpClient.enterLocalPassiveMode();
                ftpClient.setFileType(FTP.BINARY_FILE_TYPE);

                // APPROACH #1: using retrieveFile(String, OutputStream)
                String remoteFile1 = "/pub/databases/ipd/imgt/hla/alignments/"+locus+"_prot.txt";
                File downloadFile1 = new File("/data/imgtAlign/"+locus+"_prot.txt");
                OutputStream outputStream1 = new BufferedOutputStream(new FileOutputStream(downloadFile1));
                boolean success = ftpClient.retrieveFile(remoteFile1, outputStream1);
                outputStream1.close();

                if (success) {
                    System.out.println("File #1 has been downloaded successfully.");
                }else{
                    System.out.println("Something went wrong in downloading the file.");
                }

                /*// APPROACH #2: using InputStream retrieveFileStream(String)
                String remoteFile2 = "/test/song.mp3";
                File downloadFile2 = new File("D:/Downloads/song.mp3");
                OutputStream outputStream2 = new BufferedOutputStream(new FileOutputStream(downloadFile2));
                InputStream inputStream = ftpClient.retrieveFileStream(remoteFile2);
                byte[] bytesArray = new byte[4096];
                int bytesRead = -1;
                while ((bytesRead = inputStream.read(bytesArray)) != -1) {
                    outputStream2.write(bytesArray, 0, bytesRead);
                }

                success = ftpClient.completePendingCommand();
                if (success) {
                    System.out.println("File #2 has been downloaded successfully.");
                }
                outputStream2.close();
                inputStream.close();*/

            } catch (IOException ex) {
                System.out.println("Error: " + ex.getMessage());
                ex.printStackTrace();
            } finally {
                try {
                    if (ftpClient.isConnected()) {
                        ftpClient.logout();
                        ftpClient.disconnect();
                    }
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        }
    }
