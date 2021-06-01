
namespace ImageProcessingCPU
{
    partial class Form1
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.Button_LoadImage = new System.Windows.Forms.Button();
            this.MainBox = new System.Windows.Forms.PictureBox();
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            ((System.ComponentModel.ISupportInitialize)(this.MainBox)).BeginInit();
            this.SuspendLayout();
            // 
            // Button_LoadImage
            // 
            this.Button_LoadImage.Location = new System.Drawing.Point(718, 12);
            this.Button_LoadImage.Name = "Button_LoadImage";
            this.Button_LoadImage.Size = new System.Drawing.Size(75, 23);
            this.Button_LoadImage.TabIndex = 0;
            this.Button_LoadImage.Text = "New Image";
            this.Button_LoadImage.UseVisualStyleBackColor = true;
            this.Button_LoadImage.Click += new System.EventHandler(this.Button_LoadImage_Click);
            // 
            // MainBox
            // 
            this.MainBox.Location = new System.Drawing.Point(12, 12);
            this.MainBox.Name = "MainBox";
            this.MainBox.Size = new System.Drawing.Size(700, 400);
            this.MainBox.TabIndex = 1;
            this.MainBox.TabStop = false;
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(12, 415);
            this.progressBar1.Name = "progressBar1";
            this.progressBar1.Size = new System.Drawing.Size(700, 23);
            this.progressBar1.TabIndex = 2;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(7F, 15F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(799, 450);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.MainBox);
            this.Controls.Add(this.Button_LoadImage);
            this.Name = "Form1";
            this.Text = "Form1";
            this.Load += new System.EventHandler(this.Form1_Load);
            ((System.ComponentModel.ISupportInitialize)(this.MainBox)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button Button_LoadImage;
        private System.Windows.Forms.ProgressBar progressBar1;
        public System.Windows.Forms.PictureBox MainBox;
    }
}

