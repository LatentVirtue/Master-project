using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ImageProcessingCPU
{
    public partial class Form1 : Form
    {
        int canny_choice = 0;
        double lowerT = 0.1;
        double upperT = 0.3;
        public Form1()
        {
            InitializeComponent();
        }

        void SetPreview()
        {
            MainBox.Image = ImageHandler.Preview;
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void Button_LoadImage_Click(object sender, EventArgs e)
        {
            ImageHandler.Load();
            ImageHandler.Refresh();
            SetPreview();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Edge_effect(canny_choice, checkBlur.Checked);
            //MessageBox.Show(ImageHandler.Temporary.PixelFormat.ToString(), "PixelFormat: ", MessageBoxButtons.OK);
        }

        private void Button_Canny_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Canny_test(canny_choice,lowerT, upperT);
        }

        private void Button_Save_Click(object sender, EventArgs e)
        {
            if (!ImageHandler.Save())
            {
                MessageBox.Show("Image not saved.", "Mehehe", MessageBoxButtons.OK);
            }
        }

        private void canny_choice1_SelectedIndexChanged(object sender, EventArgs e)
        {
            canny_choice = canny_choice1.SelectedIndex;
        }

        private void Button_hough_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Hough();
        }

        private void label3_Click(object sender, EventArgs e)
        {

        }

        private void textBoxLowerT_TextChanged(object sender, EventArgs e)
        {
            lowerT = double.Parse(textBoxLowerT.Text);
        }

        private void textBoxUpperT_TextChanged(object sender, EventArgs e)
        {
            upperT = double.Parse(textBoxUpperT.Text);
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void textBoxNLines_TextChanged(object sender, EventArgs e)
        {
            try
            {
                ImageHandler.HoughNLines = int.Parse(textBoxNLines.Text);
            }
            catch
            {
                ImageHandler.HoughNLines = 15;
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Convolve();
        }
    }
}
