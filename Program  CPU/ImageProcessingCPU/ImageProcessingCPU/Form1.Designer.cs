
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
            this.canny_choice1 = new System.Windows.Forms.ComboBox();
            this.Button_hough = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.textBoxLowerT = new System.Windows.Forms.TextBox();
            this.textBoxUpperT = new System.Windows.Forms.TextBox();
            this.checkBlur = new System.Windows.Forms.CheckBox();
            this.textBoxNLines = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
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
            // canny_choice1
            // 
            this.canny_choice1.FormattingEnabled = true;
            this.canny_choice1.Items.AddRange(new object[] {
            "Sobel",
            "Prewitt",
            "Roberts",
            "Scharr"});
            this.canny_choice1.Location = new System.Drawing.Point(719, 92);
            this.canny_choice1.Name = "canny_choice1";
            this.canny_choice1.Size = new System.Drawing.Size(75, 23);
            this.canny_choice1.TabIndex = 6;
            this.canny_choice1.SelectedIndexChanged += new System.EventHandler(this.canny_choice1_SelectedIndexChanged);
            // 
            // Button_hough
            // 
            this.Button_hough.Location = new System.Drawing.Point(720, 179);
            this.Button_hough.Name = "Button_hough";
            this.Button_hough.Size = new System.Drawing.Size(75, 23);
            this.Button_hough.TabIndex = 7;
            this.Button_hough.Text = "Hough";
            this.Button_hough.UseVisualStyleBackColor = true;
            this.Button_hough.Click += new System.EventHandler(this.Button_hough_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(719, 71);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(54, 15);
            this.label1.TabIndex = 8;
            this.label1.Text = "Operator";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(718, 129);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(45, 15);
            this.label2.TabIndex = 9;
            this.label2.Text = "LowerT";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(718, 158);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(45, 15);
            this.label3.TabIndex = 10;
            this.label3.Text = "UpperT";
            this.label3.Click += new System.EventHandler(this.label3_Click);
            // 
            // textBoxLowerT
            // 
            this.textBoxLowerT.Location = new System.Drawing.Point(764, 121);
            this.textBoxLowerT.MaxLength = 5;
            this.textBoxLowerT.Name = "textBoxLowerT";
            this.textBoxLowerT.Size = new System.Drawing.Size(31, 23);
            this.textBoxLowerT.TabIndex = 11;
            this.textBoxLowerT.Text = "0.1";
            this.textBoxLowerT.TextChanged += new System.EventHandler(this.textBoxLowerT_TextChanged);
            // 
            // textBoxUpperT
            // 
            this.textBoxUpperT.Location = new System.Drawing.Point(764, 150);
            this.textBoxUpperT.MaxLength = 5;
            this.textBoxUpperT.Name = "textBoxUpperT";
            this.textBoxUpperT.Size = new System.Drawing.Size(31, 23);
            this.textBoxUpperT.TabIndex = 12;
            this.textBoxUpperT.Text = "0.3";
            this.textBoxUpperT.TextChanged += new System.EventHandler(this.textBoxUpperT_TextChanged);
            // 
            // checkBlur
            // 
            this.checkBlur.AutoSize = true;
            this.checkBlur.Location = new System.Drawing.Point(718, 415);
            this.checkBlur.Name = "checkBlur";
            this.checkBlur.Size = new System.Drawing.Size(47, 19);
            this.checkBlur.TabIndex = 13;
            this.checkBlur.Text = "Blur";
            this.checkBlur.UseVisualStyleBackColor = true;
            this.checkBlur.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged);
            // 
            // textBoxNLines
            // 
            this.textBoxNLines.Location = new System.Drawing.Point(718, 223);
            this.textBoxNLines.MaxLength = 5;
            this.textBoxNLines.Name = "textBoxNLines";
            this.textBoxNLines.Size = new System.Drawing.Size(31, 23);
            this.textBoxNLines.TabIndex = 14;
            this.textBoxNLines.Text = "15";
            this.textBoxNLines.TextChanged += new System.EventHandler(this.textBoxNLines_TextChanged);
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(718, 205);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(81, 15);
            this.label4.TabIndex = 15;
            this.label4.Text = "Num. of Lines";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(7F, 15F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(799, 450);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.textBoxNLines);
            this.Controls.Add(this.checkBlur);
            this.Controls.Add(this.textBoxUpperT);
            this.Controls.Add(this.textBoxLowerT);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.Button_hough);
            this.Controls.Add(this.canny_choice1);
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
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button Button_LoadImage;
        private System.Windows.Forms.ProgressBar progressBar1;
        public System.Windows.Forms.PictureBox MainBox;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button Button_Canny;
        private System.Windows.Forms.Button Button_Save;
        private System.Windows.Forms.ComboBox canny_choice1;
        private System.Windows.Forms.Button Button_hough;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox textBoxLowerT;
        private System.Windows.Forms.TextBox textBoxUpperT;
        private System.Windows.Forms.CheckBox checkBlur;
        private System.Windows.Forms.TextBox textBoxNLines;
        private System.Windows.Forms.Label label4;
    }
}

