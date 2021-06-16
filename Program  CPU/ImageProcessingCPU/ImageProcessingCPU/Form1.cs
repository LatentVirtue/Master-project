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
            MainBox.Image = ImageHandler.Edge_effect(canny_choice);
            //MessageBox.Show(ImageHandler.Temporary.PixelFormat.ToString(), "PixelFormat: ", MessageBoxButtons.OK);
        }

        private void Button_Canny_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Canny_test(canny_choice);
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
    }
}
