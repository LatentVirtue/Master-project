
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
            this.button1 = new System.Windows.Forms.Button();
            this.Button_Canny = new System.Windows.Forms.Button();
            this.Button_Save = new System.Windows.Forms.Button();
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
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(719, 388);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 3;
            this.button1.Text = "TestAlgo";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // Button_Canny
            // 
            this.Button_Canny.Location = new System.Drawing.Point(719, 42);
            this.Button_Canny.Name = "Button_Canny";
            this.Button_Canny.Size = new System.Drawing.Size(75, 23);
            this.Button_Canny.TabIndex = 4;
            this.Button_Canny.Text = "Canny";
            this.Button_Canny.UseVisualStyleBackColor = true;
            this.Button_Canny.Click += new System.EventHandler(this.Button_Canny_Click);
            // 
            // Button_Save
            // 
            this.Button_Save.Location = new System.Drawing.Point(719, 359);
            this.Button_Save.Name = "Button_Save";
            this.Button_Save.Size = new System.Drawing.Size(75, 23);
            this.Button_Save.TabIndex = 5;
            this.Button_Save.Text = "Save Image";
            this.Button_Save.UseVisualStyleBackColor = true;
            this.Button_Save.Click += new System.EventHandler(this.Button_Save_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(7F, 15F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(799, 450);
            this.Controls.Add(this.Button_Save);
            this.Controls.Add(this.Button_Canny);
            this.Controls.Add(this.button1);
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
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button Button_Canny;
        private System.Windows.Forms.Button Button_Save;
    }
}

