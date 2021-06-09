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
            MainBox.Image = ImageHandler.Rotate90();
            MessageBox.Show(ImageHandler.Temporary.PixelFormat.ToString(), "PixelFormat: ", MessageBoxButtons.OK);
        }

        private void Button_Canny_Click(object sender, EventArgs e)
        {
            MainBox.Image = ImageHandler.Canny_test();
        }

        private void Button_Save_Click(object sender, EventArgs e)
        {
            if (!ImageHandler.Save())
            {
                MessageBox.Show("Image not saved.", "Mehehe", MessageBoxButtons.OK);
            }
        }
    }
}
